const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const { Polynomial } = require("./polynomial/polynomial");
const { Evaluations } = require("./polynomial/evaluations");
const ComputeZGrandProductPolynomial = require("./grandproduct");
const readPTauHeader = require("./ptau_utils");
const { computeZHEvaluation, computeL1Evaluation } = require("./polynomial/polynomial_utils");

const logger = require("../logger.js");

module.exports = async function kzg_grandproduct_prover(pTauFilename, evalsBufferF, evalsBufferT) {
    logger.info("> KZG GRAND PRODUCT PROVER STARTED");

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve, power: nBitsPTau } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const sG1 = G1.F.n8 * 2;

    // ROUND 0. Get the settings and prepare the setup    
    evalsF = new Evaluations(evalsBufferF, curve, logger);
    evalsT = new Evaluations(evalsBufferT, curve, logger);

    // Ensure all polynomials have the same length
    if (evalsF.length() !== evalsT.length()) {
        throw new Error("Both buffers must have the same length.");
    }

    const nBits = Math.ceil(Math.log2(evalsT.length()));
    const domainSize = 2 ** nBits;

    // Ensure the polynomial has a length that is equal to domainSize
    if (evalsT.length() !== domainSize) {
        throw new Error("Polynomial length must be equal to the domain size.");
    }

    // Ensure the powers of Tau file is sufficiently large
    if (nBitsPTau < nBits) {
        throw new Error("The Powers of Tau file is not sufficiently large to commit the polynomials");
    }

    const PTau = new BigBuffer(domainSize * sG1);
    await fdPTau.readToBuffer(PTau, 0, domainSize * sG1, pTauSections[2][0].p);
    await fdPTau.close();

    logger.info("-------------------------------------");
    logger.info("  KZG GRAND PRODUCT PROVER SETTINGS");
    logger.info(`  Curve:       ${curve.name}`);
    logger.info(`  Domain size: ${domainSize}`);
    logger.info("-------------------------------------");

    let proof = {evaluations: {}, commitments: {}};
    let challenges = {};
    let polF, polT, polZ, polQ;

    const transcript = new Keccak256Transcript(curve);

    logger.info("> ROUND 1. Compute the witness polynomial commitments");
    await computeWitnessPolynomials();

    logger.info("> ROUND 2. Compute the grand-product polynomial Z ∈ 𝔽[X]");
    await ComputeZPolynomial();

    logger.info("> ROUND 3. Compute the quotient polynomial Q ∈ 𝔽[X]");
    await computeQPolynomial();

    logger.info("> ROUND 4. Compute the evaluations of the polynomials");
    computeEvaluations();

    logger.info("> ROUND 5. Compute the opening proof polynomials W𝔷, W𝔷𝛚 ∈ 𝔽[X]");
    await computeW();

    logger.info("> KZG GRAND PRODUCT PROVER FINISHED");

    return proof;

    async function computeWitnessPolynomials() {
        // Convert the evaluations to Montgomery form
        evalsF.eval = await Fr.batchToMontgomery(evalsF.eval);
        evalsT.eval = await Fr.batchToMontgomery(evalsT.eval);

        // Get the polynomials from the evaluations
        polF = await Polynomial.fromEvaluations(evalsF.eval, curve, logger);
        polT = await Polynomial.fromEvaluations(evalsT.eval, curve, logger);

        proof.commitments["F"] = await commit(polF);
        proof.commitments["T"] = await commit(polT);

        logger.info(`··· [f(x)]₁ =`, G1.toString(proof.commitments["F"]));
        logger.info(`··· [t(x)]₁ =`, G1.toString(proof.commitments["T"]));
    }

    async function ComputeZPolynomial() {
        transcript.addPolCommitment(proof.commitments["F"]);
        transcript.addPolCommitment(proof.commitments["T"]);

        challenges.gamma = transcript.getChallenge();
        logger.info("···      𝜸  =", Fr.toString(challenges.gamma));
    
        polZ = await ComputeZGrandProductPolynomial([[evalsF, evalsT]], challenges.gamma, curve, logger);
    
        proof.commitments["Z"] = await commit(polZ);
        logger.info(`··· [Z(x)]₁ =`, G1.toString(proof.commitments["Z"]));
    }

    async function computeQPolynomial() {
        transcript.addFieldElement(challenges.gamma);
        transcript.addPolCommitment(proof.commitments["Z"]);

        challenges.alpha = transcript.getChallenge();
        logger.info("···      𝜶  =", Fr.toString(challenges.alpha));

        const polZ1 = polZ.clone();
        const polL1 = await Polynomial.Lagrange1(nBits, curve, logger);
        polZ1.subScalar(Fr.one);
        await polZ1.multiply(polL1);

        const polZ21 = polZ.clone();
        await polZ21.shiftOmega();
        const polT21 = polT.clone().addScalar(challenges.gamma);
        await polZ21.multiply(polT21);

        const polZ22 = polZ.clone();
        const polF22 = polF.clone().addScalar(challenges.gamma);
        await polZ22.multiply(polF22);

        polZ21.sub(polZ22);
        polZ21.mulScalar(challenges.alpha);

        polQ = polZ1.add(polZ21);

        polQ.divZh(domainSize);

        proof.commitments["Q"] = await commit(polQ);
        logger.info(`··· [Q(x)]₁ =`, G1.toString(proof.commitments["Q"]));
    }

    function computeEvaluations() {
        transcript.addFieldElement(challenges.alpha);
        transcript.addPolCommitment(proof.commitments["Q"]);

        challenges.xi = transcript.getChallenge();
        logger.info("···      𝔷  =", Fr.toString(challenges.xi));

        proof.evaluations["fxi"] = polF.evaluate(challenges.xi);
        proof.evaluations["zxiw"] = polZ.evaluate(Fr.mul(challenges.xi, Fr.w[nBits]));

        logger.info(`···   f(𝔷)  =`, Fr.toString(proof.evaluations["fxi"]));
        logger.info(`··· Z(𝔷·𝛚)  =`, Fr.toString(proof.evaluations["zxiw"]));
    }

    async function computeW() {
        transcript.addFieldElement(challenges.xi);
        transcript.addFieldElement(proof.evaluations["fxi"]);
        transcript.addFieldElement(proof.evaluations["zxiw"]);

        challenges.v = transcript.getChallenge();
        logger.info("···      v  = ", Fr.toString(challenges.v));

        // Compute the linearisation polynomial r
        const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
        const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);

        logger.info("···  ZH(𝔷)  =", Fr.toString(ZHxi));
        logger.info("···  L₁(𝔷)  =", Fr.toString(L1xi));

        const fxi = proof.evaluations["fxi"];
        const zxiomega = proof.evaluations["zxiw"];

        const polR = polZ.clone().subScalar(Fr.one).mulScalar(L1xi);

        const polR21 = polT.clone().addScalar(challenges.gamma).mulScalar(zxiomega);
        const polR22 = polZ.clone().mulScalar(Fr.add(fxi, challenges.gamma));
        polR21.sub(polR22).mulScalar(challenges.alpha);
        polR.add(polR21);

        const polR3 = polQ.clone().mulScalar(ZHxi);
        polR.sub(polR3);

        const polWxi = polF.clone().subScalar(fxi);
        polWxi.divByXSubValue(challenges.xi);
        polWxi.mulScalar(challenges.v);

        polR.divByXSubValue(challenges.xi);
        polWxi.add(polR);

        const polWxiomega = polZ.clone().subScalar(zxiomega);
        polWxiomega.divByXSubValue(Fr.mul(challenges.xi, Fr.w[nBits]));

        proof.commitments["Wxi"] = await commit(polWxi)
        proof.commitments["Wxiw"] = await commit(polWxiomega)
        logger.info("··· [W𝔷(x)]₁   =", G1.toString(proof.commitments["Wxi"]));
        logger.info("··· [W𝔷·𝛚(x)]₁ =", G1.toString(proof.commitments["Wxiw"]));
    }

    async function commit(polynomial, name) {
        return await polynomial.multiExponentiation(PTau, name);
    }
}