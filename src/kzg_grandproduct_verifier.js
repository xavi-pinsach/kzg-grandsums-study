const { readBinFile } = require("@iden3/binfileutils");
const { Scalar } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const readPTauHeader = require("./ptau_utils");
const { computeZHEvaluation, computeL1Evaluation } = require("./polynomial/polynomial_utils");

const logger = require("../logger.js");

module.exports = async function kzg_grandproduct_verifier(proof, nBits, pTauFilename) {
    logger.info("> KZG GRAND PRODUCT VERIFIER STARTED");

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const G2 = curve.G2;

    const sG2 = G2.F.n8 * 2;
    const X2 = await fdPTau.read(sG2, pTauSections[3][0].p + sG2);
    await fdPTau.close();

    logger.info("---------------------------------------");
    logger.info("  KZG GRAND PRODUCT VERIFIER SETTINGS");
    logger.info(`  Curve:       ${curve.name}`);
    logger.info(`  Domain size: ${2 ** nBits}`);
    logger.info("---------------------------------------");

    let challenges = {};

    // STEP 1. Validate the corretness of the proof elements
    logger.info("> STEP 1. Validate [f(x)]₁,[t(x)]₁,[Z(x)]₁,[Q(x)]₁,[W𝔷(x)]₁,[W𝔷·𝛚(x)]₁ ∈ 𝔾₁");
    if(!validateCommitments()) return false;

    // STEP 2. Validate the corretness of the proof elements
    logger.info("> STEP 2. Validate f(𝔷),Z(𝔷·𝛚) ∈ 𝔽");
    if(!validateEvaluations()) return false;

    // STEP 3. Calculate challenge beta from transcript
    logger.info("> STEP 3. Compute 𝜸,𝜶,𝔷,v,u");
    computeChallenges();

    // STEP 4. Compute ZH(𝔷) and L1(𝔷)
    logger.info("> STEP 4. Compute ZH(𝔷) and L₁(𝔷)");

    const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
    const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);
    logger.info("··· ZH(𝔷) =", Fr.toString(ZHxi));
    logger.info("··· L₁(𝔷) =", Fr.toString(L1xi));

    // STEP 5. Compute r₀ := α⋅γ⋅z(𝔷·𝛚) - L₁(𝔷)
    logger.info("> STEP 5. Compute r₀");
    let r0 = Fr.mul(Fr.mul(challenges.alpha, challenges.gamma), proof.evaluations["zxiw"]);
    r0 = Fr.sub(r0, L1xi);
    logger.info("··· r₀    =", Fr.toString(r0));

    // STEP 6. Compute [D]_1 := [r'(x)]₁ + u·[Z(x)]₁, where r'(x) = r(x)-r₀
    // thus, [D]_1 = (L₁(𝔷) - α⋅(f(𝔷) + γ) + u)·[Z(x)]₁ + α⋅Z(𝔷·𝛚)·[t(x)]₁ - Z_H(𝔷)·[Q(x)]₁
    logger.info("> STEP 6. Compute [D]₁");
    let D1_12 = Fr.mul(challenges.alpha, Fr.add(proof.evaluations["fxi"], challenges.gamma));
    let D1_1 = Fr.add(Fr.sub(L1xi, D1_12), challenges.u);
    D1_1 = G1.timesFr(proof.commitments["Z"], D1_1);

    let D1_2 = Fr.mul(challenges.alpha, proof.evaluations["zxiw"]);
    D1_2 = G1.timesFr(proof.commitments["T"], D1_2);

    const D1_3 = G1.timesFr(proof.commitments["Q"], ZHxi);

    let D1 = G1.add(D1_1, G1.sub(D1_2, D1_3));
    logger.info("··· [D]₁  =", G1.toString(G1.toAffine(D1)));

    // STEP 7. Compute [F]_1 := [D]_1 + v·[f(x)]₁
    logger.info("> STEP 7. Compute [F]₁");
    let F1 = G1.timesFr(proof.commitments["F"], challenges.v);
    F1 = G1.add(D1, F1);
    logger.info("··· [F]₁  =", G1.toString(G1.toAffine(F1)));

    // STEP 8. Compute [E]_1 := (-r₀ + v·f(𝔷) + u·Z(𝔷·𝛚))·[1]_1
    logger.info("> STEP 8. Compute [E]₁");
    const E1_2 = Fr.mul(challenges.v, proof.evaluations["fxi"]);
    const E1_3 = Fr.mul(challenges.u, proof.evaluations["zxiw"]);
    let E1 = Fr.sub(Fr.add(E1_2, E1_3), r0);
    E1 = G1.timesFr(G1.one, E1);
    logger.info("··· [E]₁  =", G1.toString(G1.toAffine(E1)));

    // STEP 9. Check the pairing equation
    logger.info("> STEP 9. Check pairing equation e(-[W𝔷(x)]₁ - u·[W𝔷·𝛚(x)]₁, [x]₂)·e(𝔷·[W𝔷(x)]₁ + u𝔷ω·[W𝔷·𝛚(x)]₁ + [F]₁ - [E]₁, [1]₂) = 1");

    let A1 = G1.timesFr(proof.commitments["Wxiw"], challenges.u);
    A1 = G1.add(proof.commitments["Wxi"], A1);

    let B1 = Fr.mul(Fr.mul(challenges.u, challenges.xi), Fr.w[nBits]);
    B1 = G1.timesFr(proof.commitments["Wxiw"], B1);
    B1 = G1.add(G1.timesFr(proof.commitments["Wxi"], challenges.xi), B1);
    B1 = G1.add(B1, F1);
    B1 = G1.sub(B1, E1);
    const B2 = G2.one;

    const isValid = await curve.pairingEq(G1.neg(A1), X2, B1, B2);

    if (isValid) logger.info("> VERIFICATION OK");
    else logger.error("> VERIFICATION FAILED");

    logger.info("> KZG BASIC VERIFIER FINISHED");

    return isValid;

    function valueBelongsToField(name, value) {
        const belongs = Scalar.lt(Scalar.fromRprLE(value), Fr.p);
        if (!belongs) 
            logger.error(`··· ERROR: ${name} is not a valid field element`, Fr.toString(value));
        return belongs;
    }

    function valueBelongsToGroup1(name, value) {
        const belongs = G1.isValid(value);
        if (!belongs)
            logger.error(`··· ERROR: ${name} is not a valid G1 element`, G1.toString(value));
        return belongs;
    }

    function validateCommitments() {
        return (
            valueBelongsToGroup1("[f(x)]₁", proof.commitments["F"]) &&
            valueBelongsToGroup1("[t(x)]₁", proof.commitments["T"]) &&
            valueBelongsToGroup1("[Z(x)]₁", proof.commitments["Z"]) &&
            valueBelongsToGroup1("[Q(x)]₁", proof.commitments["Q"]) &&
            valueBelongsToGroup1("[W𝔷(x)]₁", proof.commitments["Wxi"]) &&
            valueBelongsToGroup1("[W𝔷·𝛚(x)]₁", proof.commitments["Wxiw"])
        );
    }

    function validateEvaluations() {
        return (
            valueBelongsToField("f(𝔷)", proof.evaluations["fxi"]) &&
            valueBelongsToField("Z(𝔷·𝛚)", proof.evaluations["zxiw"])
        );
    }

    function computeChallenges() {
        // STEP 1.1 Calculate challenge gamma from transcript
        const transcript = new Keccak256Transcript(curve);
        transcript.addPolCommitment(proof.commitments["F"]);
        transcript.addPolCommitment(proof.commitments["T"]);
        challenges.gamma = transcript.getChallenge();
        logger.info("··· 𝜸 = ", Fr.toString(challenges.gamma));

        // STEP 1.2 Calculate challenge alpha from transcript
        transcript.addFieldElement(challenges.gamma);
        transcript.addPolCommitment(proof.commitments["Z"]);
        challenges.alpha = transcript.getChallenge();
        logger.info("··· 𝜶 = ", Fr.toString(challenges.alpha));

        // STEP 1.3 Calculate challenge 𝔷 from transcript
        transcript.addFieldElement(challenges.alpha);
        transcript.addPolCommitment(proof.commitments["Q"]);
        challenges.xi = transcript.getChallenge();
        logger.info("··· 𝔷 = ", Fr.toString(challenges.xi));
        
        // STEP 1.4 Calculate challenge v from transcript
        transcript.addFieldElement(challenges.xi);
        transcript.addFieldElement(proof.evaluations["fxi"]);
        transcript.addFieldElement(proof.evaluations["zxiw"]);
        challenges.v = transcript.getChallenge();
        logger.info("··· v = ", Fr.toString(challenges.v));

        // STEP 1.5 Calculate challenge u from transcript
        transcript.addFieldElement(challenges.v);
        transcript.addPolCommitment(proof.commitments["Wxi"]);
        transcript.addPolCommitment(proof.commitments["Wxiw"]);
        challenges.u = transcript.getChallenge();
        logger.info("··· u = ", Fr.toString(challenges.u));
    }
};