const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const { Polynomial } = require("./polynomial/polynomial");
const { Evaluations } = require("./polynomial/evaluations");
const buildZGrandProduct = require("./grandproduct");
const readPTauHeader = require("./ptau_utils");

module.exports = async function kzg_grandproduct_prover(evalsBufferF, evalsBufferT, pTauFilename, options) {
    const logger = options.logger;

    if (logger) {
        logger.info("> KZG GRAND PRODUCT PROVER STARTED");
        logger.info("");
    }

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve, power: nBitsPTau } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;

    // STEP 0. Get the settings and prepare the setup    

    // Ensure all polynomials have the same length
    const polLen = evalsBufferF.byteLength / Fr.n8;
    if (evalsBufferF.byteLength !== evalsBufferT.byteLength) {
            throw new Error("Both buffers must have the same length.");
    }
    
    const nBits = Math.ceil(Math.log2(polLen));
    const domainSize = 2 ** nBits;

    // Ensure the polynomial has a length that is a power of two.
    if (polLen !== domainSize) {
        throw new Error("Polynomial length must be power of two.");
    }

    // Ensure the powers of Tau file is sufficiently large
    if (nBitsPTau < nBits) {
        throw new Error("Powers of Tau has not enough values for this polynomial");
    }

    const sG1 = curve.G1.F.n8 * 2;

    const PTau = new BigBuffer(domainSize * sG1);
    await fdPTau.readToBuffer(PTau, 0, domainSize * sG1, pTauSections[2][0].p);

    if (logger) {
        logger.info("-------------------------------------");
        logger.info("  KZG GRAND PRODUCT PROVER SETTINGS");
        logger.info(`  Curve:        ${curve.name}`);
        logger.info(`  #polynomials: ${polLen}`);
        logger.info("-------------------------------------");
    }

    let proof = {};
    let challenges = {};

    const transcript = new Keccak256Transcript(curve);

    logger.info("> STEP 0. Generate the polynomials from the evaluations");
    const { polF, polT } = await computePolynomials();

    logger.info("> STEP 1. Compute polynomial commitments");
    await computePolynomialCommitments();

    logger.info("> STEP 2. Compute the Z polynomial");
    const polZ = await ComputeZPolynomial();

    logger.info("> STEP 3. Compute the evaluations of the polynomials");
    computeOpenings();

    logger.info("> STEP 4. Compute Q polynomial");
    const polQ = await computeQPolynomial();

    // STEP 7. Get challenge alpha from transcript
    logger.info("> STEP 5. computeWxi");
    await computeWxi();

    if (logger) {
        logger.info("");
        logger.info("> KZG GRAND PRODUCT PROVER FINISHED");
    }

    await fdPTau.close();

    return proof;

    async function computePolynomials() {
        // Convert the evaluations to Montgomery form
        evalsBufferF = await Fr.batchToMontgomery(evalsBufferF);
        evalsBufferT = await Fr.batchToMontgomery(evalsBufferT);

        // Get the polynomials from the evaluations
        const polF = await Polynomial.fromEvaluations(evalsBufferF, curve, logger);
        const polT = await Polynomial.fromEvaluations(evalsBufferT, curve, logger);
        return { polF, polT };
    }

    async function computePolynomialCommitments() {
        proof.commitments = [];
        proof.commitments[0] = await polF.multiExponentiation(PTau, `polF`);
        proof.commitments[1] = await polT.multiExponentiation(PTau, `polT`);

        logger.info(`··· [polF(X)]_1 = `, curve.G1.toString(proof.commitments[0]));
        logger.info(`··· [polT(X)]_1 = `, curve.G1.toString(proof.commitments[1]));
    }

    async function ComputeZPolynomial() {
        transcript.addPolCommitment(proof.commitments[0]);
        transcript.addPolCommitment(proof.commitments[1]);

        challenges.beta = transcript.getChallenge();
        logger.info("··· beta = ", Fr.toString(challenges.beta));

        let polZ = [];
        polZ = await buildZGrandProduct(evalsBufferF, evalsBufferT, challenges.beta, curve, { logger });

        proof.commitmentZ = await polZ.multiExponentiation(PTau, `polZ`);
        logger.info(`··· [Z(X)]_1 = `, curve.G1.toString(proof.commitmentZ));
        return polZ;
    }

    function computeOpenings() {
        transcript.reset();
        transcript.addPolCommitment(proof.commitmentZ);

        challenges.xi = transcript.getChallenge();
        logger.info("··· xi = ", Fr.toString(challenges.xi));

        proof.evaluations = [];
        proof.evaluations[0] = polF.evaluate(challenges.xi);
        proof.evaluations[1] = polT.evaluate(challenges.xi);

        logger.info(`··· polF(xi) = `, Fr.toString(proof.evaluations[0]));
        logger.info(`··· polT(xi) = `, Fr.toString(proof.evaluations[1]));
    }

    async function computeQPolynomial() {
        transcript.reset();
        transcript.addEvaluation(proof.evaluations[0]);
        transcript.addEvaluation(proof.evaluations[1]);

        challenges.gamma = transcript.getChallenge();
        logger.info("··· gamma = ", Fr.toString(challenges.gamma));

        const evalsBufferQ = new Uint8Array(domainSize * Fr.n8);
        const evalsZ = await Evaluations.fromPolynomial(polZ, 1, curve, logger);

        let omega = Fr.one;
        for (let i = 0; i < domainSize; i++) {
            if (logger && ~i && i & (0xfff === 0)) logger.debug(`··· Q evaluation ${i}/${n}`);
            const i_n8 = i * Fr.n8;
            const i_wn8 = ((i + 1) % domainSize) * Fr.n8;

            const z_i = evalsZ.eval.slice(i_n8, i_n8 + Fr.n8);
            const z_wi = evalsZ.eval.slice(i_wn8, i_wn8 + Fr.n8);
            const f_i = evalsBufferF.slice(i_n8, i_n8 + Fr.n8);
            const t_i = evalsBufferT.slice(i_n8, i_n8 + Fr.n8);

            // IDENTITY A) L_1(x)(Z(x)-1) = 0
            const qA_i = i === 0 ? Fr.sub(z_i, Fr.one) : Fr.zero;

            //IDENTITY B) Z(Xg)(t(X)+α)−Z(X)(f(X)+α)
            const identityB1 = Fr.mul(z_wi, Fr.add(t_i, challenges.beta));
            const identityB2 = Fr.mul(z_i, Fr.add(f_i, challenges.beta));
            let qB_i = Fr.sub(identityB1, identityB2);

            //Apply gamma random factor
            qB_i = Fr.mul(qB_i, challenges.gamma);

            const q_i = Fr.add(qA_i, qB_i);

            evalsBufferQ.set(q_i, i_n8);

            omega = Fr.mul(omega, Fr.w[nBits]);
        }

        if (logger) logger.debug("··· Interpolating Q polynomial");
        const polQ = await Polynomial.fromEvaluations(evalsBufferQ, curve, logger);
        polQ.divZh();

        proof.commitmentQ = await polQ.multiExponentiation(PTau, `polQ`);
        logger.info(`··· [Q(X)]_1 = `, curve.G1.toString(proof.commitmentQ));
        return polQ;
    }

    async function computeWxi() {
        transcript.reset();
        transcript.addEvaluation(challenges.gamma);

        challenges.alpha = transcript.getChallenge();
        logger.info("··· alpha = ", Fr.toString(challenges.alpha));

        let polWxi = new Polynomial(new Uint8Array(polLen * Fr.n8), curve, logger);

        polQ.divByXSubValue(challenges.xi);
        polWxi.add(polQ);

        let currentAlpha = challenges.alpha;
        const pols = [polF, polT];
        for (let i = 0; i < pols.length; i++) {
            pols[i].subScalar(proof.evaluations[i]);
            pols[i].divByXSubValue(challenges.xi);
            pols[i].mulScalar(currentAlpha);

            polWxi.add(pols[i]);
            currentAlpha = Fr.mul(currentAlpha, challenges.alpha);
        }

        proof.commitWxi = await polWxi.multiExponentiation(PTau, "Wxi");
        logger.info("··· [Wxi(X)]_1 = ", curve.G1.toString(proof.commitWxi));
    }
}