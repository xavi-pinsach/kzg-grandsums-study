const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const { Polynomial } = require("./polynomial/polynomial");
const { Evaluations } = require("./polynomial/evaluations");
const ComputeSGrandSumPolynomial = require("./grandsum");
const readPTauHeader = require("./ptau_utils");
const { computeZHEvaluation, computeL1Evaluation } = require("./polynomial/polynomial_utils");

const logger = require("../logger.js");

module.exports = async function kzg_grandsum_prover(evalsBufferF, evalsBufferT, pTauFilename) {
    logger.info("> KZG GRAND SUM PROVER STARTED");

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

    const PTau = new BigBuffer(domainSize * 2 * sG1);
    await fdPTau.readToBuffer(PTau, 0, domainSize * 2 * sG1, pTauSections[2][0].p);
    await fdPTau.close();

    logger.info("-------------------------------------");
    logger.info("  KZG GRAND SUM PROVER SETTINGS");
    logger.info(`  Curve:       ${curve.name}`);
    logger.info(`  Domain size: ${domainSize}`);
    logger.info("-------------------------------------");

    let proof = {evaluations: {}, commitments: {}};
    let challenges = {};
    let polF, polT, polS, polQ;

    const transcript = new Keccak256Transcript(curve);

    logger.info("> ROUND 1. Generate the witness polynomials f,t ∈ 𝔽[X] from the evaluations");
    await computeWitnessPolynomials();

    logger.info("> ROUND 2. Compute the grand-sum polynomial S ∈ 𝔽[X]");
    await ComputeSPolynomial();

    logger.info("> ROUND 3. Compute the quotient polynomial Q ∈ 𝔽[X]");
    await computeQPolynomial();

    logger.info("> ROUND 4. Compute the evaluations of the polynomials");
    computeEvaluations();

    logger.info("> ROUND 5. Compute the opening proof polynomials W𝔷, W𝔷𝛚 ∈ 𝔽[X]");
    await computeW();

    logger.info("");
    logger.info("> KZG GRAND SUM PROVER FINISHED");

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

    async function ComputeSPolynomial() {
        transcript.addPolCommitment(proof.commitments["F"]);
        transcript.addPolCommitment(proof.commitments["T"]);

        challenges.gamma = transcript.getChallenge();
        logger.info("···      𝜸  =", Fr.toString(challenges.gamma));

        polS = await ComputeSGrandSumPolynomial([[evalsF, evalsT]], challenges.gamma, curve);

        proof.commitments["S"] = await commit(polS);
        logger.info(`··· [S(x)]₁ =`, G1.toString(proof.commitments["S"]));
    }

    async function computeQPolynomial() {
        transcript.addFieldElement(challenges.gamma);
        transcript.addPolCommitment(proof.commitments["S"]);

        challenges.alpha = transcript.getChallenge();
        logger.info("···      𝜶  =", Fr.toString(challenges.alpha));

        const polS1 = polS.clone();
        const polL1 = await Polynomial.Lagrange1(nBits, curve, logger);
        await polS1.multiply(polL1);

        const polS21 = polS.clone();
        await polS21.shiftOmega();
        polS21.sub(polS);

        const polS22 = polF.clone().addScalar(challenges.gamma);
        const polS23 = polT.clone().addScalar(challenges.gamma);

        await polS21.multiply(polS22);
        await polS21.multiply(polS23);

        polS21.add(polF);
        polS21.sub(polT);

        await polS21.mulScalar(challenges.alpha);

        polQ = polS1.add(polS21);

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
        proof.evaluations["txi"] = polT.evaluate(challenges.xi);
        proof.evaluations["sxiw"] = polS.evaluate(Fr.mul(challenges.xi, Fr.w[nBits]));

        logger.info(`···   f(𝔷)  =`, Fr.toString(proof.evaluations["fxi"]));
        logger.info(`···   t(𝔷)  =`, Fr.toString(proof.evaluations["txi"]));
        logger.info(`··· S(𝔷·𝛚)  =`, Fr.toString(proof.evaluations["sxiw"]));
    }

    async function computeW() {
        transcript.addFieldElement(challenges.xi);
        transcript.addFieldElement(proof.evaluations["fxi"]);
        transcript.addFieldElement(proof.evaluations["txi"]);
        transcript.addFieldElement(proof.evaluations["sxiw"]);

        challenges.v = transcript.getChallenge();
        logger.info("···      v  = ", Fr.toString(challenges.v));

        const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
        const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);

        logger.info("···  ZH(𝔷)  =", Fr.toString(ZHxi));
        logger.info("···  L₁(𝔷)  =", Fr.toString(L1xi));

        const fxi = proof.evaluations["fxi"];
        const txi = proof.evaluations["txi"];
        const sxiomega = proof.evaluations["sxiw"];

        const polR = polS.clone().mulScalar(L1xi);

        const polR2 = polS.clone().mulScalar(Fr.negone).addScalar(sxiomega);
        polR2.mulScalar(Fr.add(fxi, challenges.gamma));
        polR2.mulScalar(Fr.add(txi, challenges.gamma));
        polR2.addScalar(Fr.sub(fxi, txi));
        polR2.mulScalar(challenges.alpha);

        polR.add(polR2);

        const polR3 = polQ.clone().mulScalar(ZHxi);
        polR.sub(polR3);

        // The following can be optimized by using Homer's rule
        const polWxi = polF.clone().subScalar(fxi);
        polWxi.mulScalar(challenges.v);
        polWxi.divByXSubValue(challenges.xi);

        const polWxi2 = polT.clone().subScalar(txi);
        polWxi2.mulScalar(Fr.square(challenges.v));
        polWxi2.divByXSubValue(challenges.xi);
        polWxi.add(polWxi2);

        polR.divByXSubValue(challenges.xi);
        polWxi.add(polR);

        const polWxiomega = polS.clone().subScalar(sxiomega);
        polWxiomega.divByXSubValue(Fr.mul(challenges.xi, Fr.w[nBits]));

        proof.commitments["Wxi"] = await commit(polWxi);
        proof.commitments["Wxiw"] = await commit(polWxiomega);
        logger.info("··· [W𝔷(x)]₁   =", G1.toString(proof.commitments["Wxi"]));
        logger.info("··· [W𝔷·𝛚(x)]₁ =", G1.toString(proof.commitments["Wxiw"]));
    }

    async function commit(polynomial, name) {
        return await polynomial.multiExponentiation(PTau, name);
    }
}