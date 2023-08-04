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

    const PTau = new BigBuffer(domainSize * sG1);
    await fdPTau.readToBuffer(PTau, 0, domainSize * sG1, pTauSections[2][0].p);

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

        let extension = 4;
        const extensionBits = Math.log2(extension);

        const evalsQ = new Evaluations(new Uint8Array(domainSize * extension * Fr.n8), curve, logger);
        const evalsS = await Evaluations.fromPolynomial(polS, extension, curve, logger);
        const evalsFExt = await Evaluations.fromPolynomial(polF, extension, curve, logger);
        const evalsTExt = await Evaluations.fromPolynomial(polT, extension, curve, logger);

        let omega = Fr.one;
        for (let i = 0; i < domainSize * extension; i++) {
            if (logger && ~i && i & (0xfff === 0)) logger.debug(`··· Q evaluation ${i}/${n}`);

            const s_i = evalsS.getEvaluation(i);
            const s_wi = evalsS.getEvaluation((i + extension) % (domainSize * extension));
            const f_i = evalsFExt.getEvaluation(i);
            const t_i = evalsTExt.getEvaluation(i);

            const ZH_i = computeZHEvaluation(curve, omega, nBits);
            const L1_i = computeL1Evaluation(curve, omega, ZH_i, nBits);

            // IDENTITY A) L_1(x)S(x) = 0
            const qA_i = Fr.mul(L1_i, s_i);

            //IDENTITY B) (S(Xg) - S(X))·(f(X) + 𝜸)·(t(X) + 𝜸) + f(X) - t(X)
            const qB11_i = Fr.sub(s_wi, s_i);
            const qB12_i = Fr.add(f_i, challenges.gamma);
            const qB13_i = Fr.add(t_i, challenges.gamma);
            const qB1_i = Fr.mul(qB11_i, Fr.mul(qB12_i, qB13_i));

            let qB_i = Fr.add(qB1_i, Fr.sub(f_i, t_i));
            //Apply alpha random factor
            qB_i = Fr.mul(qB_i, challenges.alpha);

            const q_i = Fr.add(qA_i, qB_i);
            evalsQ.setEvaluation(i, q_i);

            // Compute next omega
            omega = Fr.mul(omega, Fr.w[nBits + extensionBits]);
        }

        if (logger) logger.debug("··· Interpolating Q polynomial");
        const polQ = await Polynomial.fromEvaluations(evalsQ.eval, curve, logger);
        polQ.divZh(domainSize, extension);

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

        // Compute the linearisation polynomial r
        const evalsR = new Evaluations(new Uint8Array(domainSize * Fr.n8), curve, logger);
        const evalsS = await Evaluations.fromPolynomial(polS, 1, curve, logger);
        const evalsQ = await Evaluations.fromPolynomial(polQ, 1, curve, logger);

        const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
        const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);

        logger.info("···  ZH(𝔷)  =", Fr.toString(ZHxi));
        logger.info("···  L₁(𝔷)  =", Fr.toString(L1xi));

        const fxi = proof.evaluations["fxi"];
        const txi = proof.evaluations["txi"];
        const sxiomega = proof.evaluations["sxiw"];

        let omega = Fr.one; // I think this lookp sould be until 2n at least
        for (let i = 0; i < domainSize; i++) {
            if (logger && ~i && i & (0xfff === 0)) logger.debug(`··· R evaluation ${i}/${n}`);

            const s_i = evalsS.getEvaluation(i);
            const q_i = evalsQ.getEvaluation(i);

            // Factor 1. S(x)L_1(𝔷)
            const rA_i = Fr.mul(s_i, L1xi);

            // Factor 2. (S(𝔷·𝛚) - S(X))·(f(𝔷) + 𝜸)·(t(𝔷) + 𝜸) + f(𝔷) - t(𝔷)
            const rB11_i = Fr.sub(sxiomega, s_i);
            const rB12_i = Fr.add(fxi, challenges.gamma);
            const rB13_i = Fr.add(txi, challenges.gamma);
            const rB1_i = Fr.mul(rB11_i, Fr.mul(rB12_i, rB13_i));

            let rB_i = Fr.add(rB1_i, Fr.sub(fxi, txi));
            rB_i = Fr.mul(rB_i, challenges.alpha);

            // Factor 3. ZH(𝔷)·q(X)
            const rC_i = Fr.mul(ZHxi, q_i);

            let r_i = Fr.add(rA_i, rB_i);
            r_i = Fr.sub(r_i, rC_i);

            evalsR.setEvaluation(i, r_i);

            omega = Fr.mul(omega, Fr.w[nBits]);
        }

        if (logger) logger.debug("··· Interpolating r polynomial");
        const polR = await Polynomial.fromEvaluations(evalsR.eval, curve, logger);

        challenges.v = transcript.getChallenge();
        logger.info("···      v  = ", Fr.toString(challenges.v));

        let polWxi = new Polynomial(new Uint8Array(domainSize * Fr.n8), curve, logger);
        let polWxiomega = new Polynomial(new Uint8Array(domainSize * Fr.n8), curve, logger);


        polF.subScalar(fxi);
        polF.mulScalar(challenges.v);
        polF.divByXSubValue(challenges.xi);
        polWxi.add(polF);

        polT.subScalar(txi);
        polT.mulScalar(Fr.square(challenges.v));
        polT.divByXSubValue(challenges.xi);
        polWxi.add(polF);

        polR.divByXSubValue(challenges.xi);
        polWxi.add(polR);

        polS.subScalar(sxiomega);
        polWxiomega.add(polS);
        polWxiomega.divByXSubValue(Fr.mul(challenges.xi, Fr.w[nBits]));

        proof.commitments["Wxi"] = await polWxi.multiExponentiation(PTau, "Wxi");
        proof.commitments["Wxiw"] = await polWxiomega.multiExponentiation(PTau, "Wxiomega");
        logger.info("··· [W𝔷(x)]₁   =", G1.toString(proof.commitments["Wxi"]));
        logger.info("··· [W𝔷·𝛚(x)]₁ =", G1.toString(proof.commitments["Wxiw"]));
    }
}