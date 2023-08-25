const { readBinFile } = require("@iden3/binfileutils");
const { Scalar } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const readPTauHeader = require("./ptau_utils");
const { computeZHEvaluation, computeL1Evaluation } = require("./polynomial/polynomial_utils");

const logger = require("../logger.js");

module.exports = async function mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits, nPols = 1) {
    logger.info("> KZG GRAND SUM VERIFIER STARTED");

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const G2 = curve.G2;

    const sG2 = G2.F.n8 * 2;
    const X2 = await fdPTau.read(sG2, pTauSections[3][0].p + sG2);
    await fdPTau.close();

    logger.info("---------------------------------------");
    logger.info("  KZG GRAND SUM VERIFIER SETTINGS");
    logger.info(`  Curve:       ${curve.name}`);
    logger.info(`  Domain size: ${2 ** nBits}`);
    logger.info("---------------------------------------");

    let challenges = {};
    const isVector = nPols > 1;
    let step = 1;

    let pols = "";
    if (isVector) {
        for (let i = 1; i <= nPols; i++) {
            pols += `[f${i}(x)]₁,[t${i}(x)]₁,`;
        }
    }
    logger.info(`> STEP ${step}. Validate ${pols}[f(x)]₁,[t(x)]₁,[S(x)]₁,[Q(x)]₁,[W𝔷(x)]₁,[W𝔷·𝛚(x)]₁ ∈ 𝔾₁`);

    if(!validateCommitments()) return false;
    ++step;

    logger.info(`> STEP ${step}. Validate f(𝔷),t(𝔷),S(𝔷·𝛚) ∈ 𝔽`);
    if(!validateEvaluations()) return false;
    ++step;

    if (isVector) {
        logger.info(`> STEP ${step}. Compute 𝛃,𝜸,𝜶,𝔷,v,u`);
    } else {
        logger.info(`> STEP ${step}. Compute 𝜸,𝜶,𝔷,v,u`);
    }
    computeChallenges();
    ++step;

    if (isVector) {
        logger.info(`> STEP ${step}. Check that the fi and ti are correctly related to the f and t polynomials, respectively`);
        if(!validatePolynomialRelation()) return false;
        ++step;
    }

    logger.info(`> STEP ${step}. Compute ZH(𝔷) and L₁(𝔷)`);
    const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
    const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);
    logger.info("··· ZH(𝔷) =", Fr.toString(ZHxi));
    logger.info("··· L₁(𝔷) =", Fr.toString(L1xi));
    ++step;

    logger.info(`> STEP ${step}. Compute r₀ = α⋅[ S(𝔷·𝛚)·(f(𝔷) + γ)·(t(𝔷) + γ) + f(𝔷) - t(𝔷) ]`);
    const r0_2 = Fr.add(proof.evaluations["fxi"], challenges.gamma);
    const r0_3 = Fr.add(proof.evaluations["txi"], challenges.gamma);
    let r0 = Fr.mul(proof.evaluations["sxiw"], Fr.mul(r0_2, r0_3));
    r0 = Fr.add(r0, Fr.sub(proof.evaluations["fxi"], proof.evaluations["txi"]))
    r0 = Fr.mul(r0, challenges.alpha);
    logger.info("··· r₀    =", Fr.toString(r0));
    ++step;

    logger.info(`> STEP ${step}. Compute [D]_1 = (L₁(𝔷) - α⋅(f(𝔷) + γ)·(t(𝔷) + γ) + u)·[S(x)]₁ - Z_H(𝔷)·[Q(x)]₁`);
    let D1_12 = Fr.mul(challenges.alpha, Fr.add(proof.evaluations["fxi"], challenges.gamma));
    D1_12 = Fr.mul(D1_12, Fr.add(proof.evaluations["txi"], challenges.gamma));
    let D1_1 = Fr.add(Fr.sub(L1xi, D1_12), challenges.u);
    D1_1 = G1.timesFr(proof.commitments["S"], D1_1);

    const D1_2 = G1.timesFr(proof.commitments["Q"], ZHxi);

    let D1 = G1.sub(D1_1,D1_2);
    logger.info("··· [D]₁  =", G1.toString(G1.toAffine(D1)));
    ++step;

    logger.info(`> STEP ${step}. Compute [F]₁ = [D]_1 + v·[f(x)]₁ + v^2·[t(x)]₁`);
    let F1 = G1.timesFr(proof.commitments["F"], challenges.v);
    F1 = G1.add(F1, G1.timesFr(proof.commitments["T"], Fr.square(challenges.v)));
    F1 = G1.add(F1, D1);
    logger.info("··· [F]₁  =", G1.toString(G1.toAffine(F1)));
    ++step;

    logger.info(`> STEP ${step}. Compute [E]₁ = (-r₀ + v·f(𝔷) + v^2·t(𝔷) + u·S(𝔷·𝛚))·[1]_1`);
    const E1_2 = Fr.mul(challenges.v, proof.evaluations["fxi"]);
    const E1_3 = Fr.mul(Fr.square(challenges.v), proof.evaluations["txi"]);
    const E1_4 = Fr.mul(challenges.u, proof.evaluations["sxiw"]);
    let E1 = Fr.sub(Fr.add(E1_2, Fr.add(E1_3, E1_4)), r0);
    E1 = G1.timesFr(G1.one, E1);
    logger.info("··· [E]₁  =", G1.toString(G1.toAffine(E1)));
    ++step;

    logger.info(`> STEP ${step}. Check pairing equation e(-[W𝔷(x)]₁ - u·[W𝔷·𝛚(x)]₁, [x]₂)·e(𝔷·[W𝔷(x)]₁ + u𝔷ω·[W𝔷·𝛚(x)]₁ + [F]₁ - [E]₁, [1]₂) = 1`);
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

    logger.info("> KZG GRAND SUM VERIFIER FINISHED");

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
        if (isVector) {
            for (let i = 0; i < nPols; i++) {
                if (!valueBelongsToGroup1(`[f${i}(x)]₁`, proof.commitments[`F${i}`])) return false;
                if (!valueBelongsToGroup1(`[t${i}(x)]₁`, proof.commitments[`T${i}`])) return false;
            }
        }

        return (
            valueBelongsToGroup1("[f(x)]₁", proof.commitments["F"]) &&
            valueBelongsToGroup1("[t(x)]₁", proof.commitments["T"]) &&
            valueBelongsToGroup1("[S(x)]₁", proof.commitments["S"]) &&
            valueBelongsToGroup1("[Q(x)]₁", proof.commitments["Q"]) &&
            valueBelongsToGroup1("[W𝔷(x)]₁", proof.commitments["Wxi"]) &&
            valueBelongsToGroup1("[W𝔷·𝛚(x)]₁", proof.commitments["Wxiw"])
        );
    }

    function validateEvaluations() {
        return (
            valueBelongsToField("f(𝔷)", proof.evaluations["fxi"]) &&
            valueBelongsToField("t(𝔷)", proof.evaluations["txi"]) &&
            valueBelongsToField("S(𝔷·𝛚)", proof.evaluations["sxiw"])
        );
    }

    function computeChallenges() {
        const transcript = new Keccak256Transcript(curve);

        if (isVector) {
            // STEP 1.1 Calculate challenge beta from transcript
            for (let i = 0; i < nPols; i++) {
                transcript.addPolCommitment(proof.commitments[`F${i}`]);
                transcript.addPolCommitment(proof.commitments[`T${i}`]);
            }
            challenges.beta = transcript.getChallenge();
            logger.info("··· 𝛽 =", Fr.toString(challenges.beta));
        }
        
        // STEP 1.2 Calculate challenge gamma from transcript
        transcript.addPolCommitment(proof.commitments["F"]);
        transcript.addPolCommitment(proof.commitments["T"]);
        challenges.gamma = transcript.getChallenge();
        logger.info("··· 𝜸 =", Fr.toString(challenges.gamma));

        // STEP 1.3 Calculate challenge alpha from transcript
        transcript.addFieldElement(challenges.gamma);
        transcript.addPolCommitment(proof.commitments["S"]);
        challenges.alpha = transcript.getChallenge();
        logger.info("··· 𝜶 =", Fr.toString(challenges.alpha));

        // STEP 1.4 Calculate challenge 𝔷 from transcript
        transcript.addFieldElement(challenges.alpha);
        transcript.addPolCommitment(proof.commitments["Q"]);
        challenges.xi = transcript.getChallenge();
        logger.info("··· 𝔷 =", Fr.toString(challenges.xi));
        
        // STEP 1.5 Calculate challenge v from transcript
        transcript.addFieldElement(challenges.xi);
        transcript.addFieldElement(proof.evaluations["fxi"]);
        transcript.addFieldElement(proof.evaluations["txi"]);
        transcript.addFieldElement(proof.evaluations["sxiw"]);
        challenges.v = transcript.getChallenge();
        logger.info("··· v =", Fr.toString(challenges.v));

        // STEP 1.6 Calculate challenge u from transcript
        transcript.addFieldElement(challenges.v);
        transcript.addPolCommitment(proof.commitments["Wxi"]);
        transcript.addPolCommitment(proof.commitments["Wxiw"]);
        challenges.u = transcript.getChallenge();
        logger.info("··· u =", Fr.toString(challenges.u));
    }

    function validatePolynomialRelation() {
        firstPol = proof.commitments[`F${nPols - 1}`];
        secondPol = proof.commitments[`T${nPols - 1}`];
        for (let i = nPols - 2; i >= 0; i--) {
            firstPol = G1.add(G1.timesFr(firstPol, challenges.beta), proof.commitments[`F${i}`]);
            secondPol = G1.add(G1.timesFr(secondPol, challenges.beta), proof.commitments[`T${i}`]);
        }

        return G1.eq(firstPol, proof.commitments["F"]) && G1.eq(secondPol, proof.commitments["T"]);
    }
};