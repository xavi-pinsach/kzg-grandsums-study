const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const { Polynomial } = require("./polynomial/polynomial");
const { Evaluations } = require("./polynomial/evaluations");
const buildSGrandsum = require("./grandsum");
const readPTauHeader = require("./ptau_utils");
const { computeZHEvaluation, computeL1Evaluation } = require("./polynomial/polynomial_utils");

module.exports = async function kzg_grandsum_prover(evalsBufferF, evalsBufferT, pTauFilename, options) {
    const logger = options.logger;

    if (logger) {
        logger.info("> KZG GRAND SUM PROVER STARTED");
        logger.info("");
    }

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve, power: nBitsPTau } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const sG1 = G1.F.n8 * 2;

    // STEP 0. Get the settings and prepare the setup    
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
        throw new Error("Powers of Tau has not enough values for this polynomial");
    }

    const PTau = new BigBuffer(domainSize * 2 * sG1);
    await fdPTau.readToBuffer(PTau, 0, domainSize * 2 * sG1, pTauSections[2][0].p);

    if (logger) {
        logger.info("-------------------------------------");
        logger.info("  KZG GRAND SUM PROVER SETTINGS");
        logger.info(`  Curve:       ${curve.name}`);
        logger.info(`  Domain size: ${domainSize}`);
        logger.info("-------------------------------------");
    }

    let proof = {evaluations: {}, commitments: {}};
    let challenges = {};

    const transcript = new Keccak256Transcript(curve);

    logger.info("> STEP 0. Generate the witness polynomials f,t ∈ 𝔽[X] from the evaluations");
    const { polF, polT } = await computeWitnessPolynomials();

    logger.info("> STEP 1. Compute the witness polynomial commitments");
    await computeWitnessPolsCommitments();

    logger.info("> STEP 2. Compute the grand-sum polynomial S ∈ 𝔽[X]");
    const polS = await ComputeSPolynomial();

    logger.info("> STEP 3. Compute the quotient polynomial Q ∈ 𝔽[X]");
    const polQ = await computeQPolynomial();

    logger.info("> STEP 4. Compute the evaluations of the polynomials");
    computeEvaluations();

    logger.info("> STEP 5. Compute the opening proof polynomials W𝔷, W𝔷𝛚 ∈ 𝔽[X]");
    await computeW();

    if (logger) {
        logger.info("");
        logger.info("> KZG GRAND SUM PROVER FINISHED");
    }

    await fdPTau.close();

    return proof;

    async function computeWitnessPolynomials() {
        // Convert the evaluations to Montgomery form
        evalsF.eval = await Fr.batchToMontgomery(evalsF.eval);
        evalsT.eval = await Fr.batchToMontgomery(evalsT.eval);

        // Get the polynomials from the evaluations
        const polF = await Polynomial.fromEvaluations(evalsF.eval, curve, logger);
        const polT = await Polynomial.fromEvaluations(evalsT.eval, curve, logger);
        return { polF, polT };
    }

    async function computeWitnessPolsCommitments() {
        proof.commitments["F"] = await polF.multiExponentiation(PTau, `polF`);
        proof.commitments["T"] = await polT.multiExponentiation(PTau, `polT`);

        logger.info(`··· [f(x)]₁ =`, G1.toString(proof.commitments["F"]));
        logger.info(`··· [t(x)]₁ =`, G1.toString(proof.commitments["T"]));
    }

    async function ComputeSPolynomial() {
        transcript.addPolCommitment(proof.commitments["F"]);
        transcript.addPolCommitment(proof.commitments["T"]);

        challenges.gamma = transcript.getChallenge();
        logger.info("···      𝜸  =", Fr.toString(challenges.gamma));

        let polS = [];
        polS = await buildSGrandsum(evalsF, evalsT, challenges.gamma, curve, { logger });

        proof.commitments["S"] = await polS.multiExponentiation(PTau, `polS`);
        logger.info(`··· [S(x)]₁ =`, G1.toString(proof.commitments["S"]));
        return polS;
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

        const polQ = polS1.add(polS21);

        polQ.divZh(domainSize);

        proof.commitments["Q"] = await polQ.multiExponentiation(PTau, `polQ`);
        logger.info(`··· [Q(x)]₁ =`, G1.toString(proof.commitments["Q"]));
        return polQ;
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

        proof.commitments["Wxi"] = await polWxi.multiExponentiation(PTau, "Wxi");
        proof.commitments["Wxiw"] = await polWxiomega.multiExponentiation(PTau, "Wxiomega");
        logger.info("··· [W𝔷(x)]₁   =", G1.toString(proof.commitments["Wxi"]));
        logger.info("··· [W𝔷·𝛚(x)]₁ =", G1.toString(proof.commitments["Wxiw"]));
    }
}