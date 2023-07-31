const { readBinFile } = require("@iden3/binfileutils");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const readPTauHeader = require("./ptau_utils");

module.exports = async function kzg_grandproduct_verifier(proof, pTauFilename, options) {
    const logger = options.logger;

    if (logger) {
        logger.info("> KZG BASIC VERIFIER STARTED");
        logger.info("");
    }

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve } = await readPTauHeader(fdPTau, pTauSections);

    const nPols = proof.commitments.length;
    if (logger) {
        logger.info("---------------------------------------");
        logger.info("  KZG GRAND PRODUCT VERIFIER SETTINGS");
        logger.info(`  Curve:        ${curve.name}`);
        logger.info(`  #polynomials: ${nPols}`);
        logger.info("---------------------------------------");
    }

    let challenges = {};

    // STEP 1. Calculate challenge beta from transcript
    logger.info("> STEP 1.1 Compute challenge beta");
    const transcript = new Keccak256Transcript(curve);
    for (let i = 0; i < nPols; i++) {
        transcript.addPolCommitment(proof.commitments[i]);
    }
    challenges.beta = transcript.getChallenge();
    logger.info("··· beta = ", curve.Fr.toString(challenges.beta));

    // STEP 1. Calculate challenge xi from transcript
    logger.info("> STEP 1.2 Compute challenge xi");
    transcript.reset();
    for (let i = 0; i < nPols; i++) {
        transcript.addPolCommitment(proof.commitments[i]);
    }
    challenges.xi = transcript.getChallenge();
    logger.info("··· xi = ", curve.Fr.toString(challenges.xi));

    // STEP 2. Calculate challenge alpha from transcript
    logger.info("> STEP 2. Compute challenge alpha");
    transcript.reset();
    for (let i = 0; i < nPols; i++) {
        transcript.addEvaluation(proof.evaluations[i]);
    }
    challenges.alpha = transcript.getChallenge();
    logger.info("··· alpha = ", curve.Fr.toString(challenges.alpha));
    


    // STEP 3. Compute r0
    const r0 = 

    // STEP 3. Compute [F]_1
    // let currentAlpha = curve.Fr.one;
    let F = proof.commitmentQ;
    // for(let i = 0; i < nPols; i++) {
    //     F = curve.G1.add(F, curve.G1.timesFr(proof.commitments[i], currentAlpha));
    //     currentAlpha = curve.Fr.mul(currentAlpha, challenges.alpha);
    // }
    F = curve.G1.add(F, curve.G1.timesFr(proof.commitments[0], challenges.alpha));

    // STEP 4. Compute [E]_1
    currentAlpha = curve.Fr.one;
    let e = curve.Fr.zero;
    for(let i = 0; i < nPols; i++) {
        e = curve.Fr.add(e, curve.Fr.mul(proof.evaluations[i], currentAlpha));
        currentAlpha = curve.Fr.mul(currentAlpha, challenges.alpha);
    }
    const E = curve.G1.timesFr(curve.G1.one, e);

    // STEP 3. Check the pairing equation
    logger.info("> STEP 3. Check pairing");

    const A1 = proof.commitWxi

    const sG2 = curve.G2.F.n8 * 2;
    const S_2 = await fdPTau.read(sG2, pTauSections[3][0].p + sG2);
    const xi_2 = curve.G2.timesFr(curve.G2.one, challenges.xi);
    const A2 = curve.G2.sub(S_2, xi_2);

    const B1 = curve.G1.sub(F, E);
    const B2 = curve.G2.one;

    const isValid = await curve.pairingEq(A1, A2, B1, B2);

    if (logger) {
        if (isValid) {
            logger.info("> VERIFICATION OK");
        } else {
            logger.error("> VERIFICATION FAILED");
        }
    }

    if (logger) {
        logger.info("");
        logger.info("> KZG BASIC VERIFIER FINISHED");
    }

    await fdPTau.close();

    return isValid;
};

// TODO: Why some are put in challenges and others are returned?
function calculateL1andZHEvaluation(curve, challenges, vk) {
    const Fr = curve.Fr;

    let xin = challenges.xi;
    let domainSize = 1;
    for (let i=0; i<vk.power; i++) {
        xin = Fr.square(xin);
        domainSize *= 2;
    }
    // challenges.xin = xin;

    // challenges.zh = Fr.sub(xin, Fr.one);
    const ZHxi = Fr.sub(xin, Fr.one);

    // const L = [];

    // const n = Fr.e(domainSize);
    // let w = Fr.one;
    // for (let i=0; i<=Math.max(1, vk.nPublic); i++) {
    //     L[i] = Fr.div(Fr.mul(w, challenges.zh), Fr.mul(n, Fr.sub(challenges.xi, w)));
    //     w = Fr.mul(w, Fr.w[vk.power]);
    // }

    // return L;

    const n = Fr.e(domainSize);
    const w = Fr.w[vk.power];
    const L1xi = Fr.div(Fr.mul(w, ZHxi), Fr.mul(n, Fr.sub(challenges.xi, w)));

    return { ZHxi, L1xi };
}
