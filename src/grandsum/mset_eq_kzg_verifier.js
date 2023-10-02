const { readBinFile } = require("@iden3/binfileutils");
const { Scalar } = require("ffjavascript");
const { Keccak256Transcript } = require("../Keccak256Transcript");
const { computeZHEvaluation, computeL1Evaluation } = require("../polynomial/polynomial_utils");
const { readPTauHeader } = require("../ptau_utils");

const logger = require("../../logger.js");

module.exports = async function mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits) {
    logger.info("> MULTISET EQUALITY KZG GRAND-SUM VERIFIER STARTED");
    
    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const G2 = curve.G2;

    const sG2 = G2.F.n8 * 2;
    const X2 = await fdPTau.read(sG2, pTauSections[3][0].p + sG2);
    await fdPTau.close();

    // Obtain the number of polynomials subject to the argument (proof.commitments.Fi, where i ∈ [n]) from the proof.
    const nFiCommitments = Object.keys(proof.commitments).filter(k => k.match(/^F\d/)).length;
    const nPols = nFiCommitments > 0 ? nFiCommitments : 1;
    const isVector = nPols > 1;

    // Obtain the selectiveness of the argument from the proof
    const isSelected = Object.keys(proof.commitments).filter(k => k.match(/^selF/)).length === 1;

    logger.info("---------------------------------------");
    logger.info("  MULTISET EQUALITY KZG GRAND-SUM VERIFIER SETTINGS");
    logger.info(`  Curve:        ${curve.name}`);
    logger.info(`  Domain size:  ${2 ** nBits}`);
    logger.info(`  #polynomials: ${nPols}`);
    logger.info(`  Selectors:    ${isSelected ? "Yes" : "No"}`);
    logger.info("---------------------------------------");

    let challenges = {};
    let step = 1;

    let pols = "";
    if (isVector) {
        for (let i = 0; i < nPols; i++) {
            pols += `[f${i+1}(x)]₁,[t${i+1}(x)]₁,`;
        }
    }
    if (isSelected) pols += "[fsel(x)]₁,[tsel(x)]₁,";
    logger.info(`> STEP ${step}. Validate ${pols}[S(x)]₁,[Q(x)]₁,[W𝔷(x)]₁,[W𝔷·𝛚(x)]₁ ∈ 𝔾₁`);

    if(!validateCommitments()) return false;
    ++step;

    let evals = "";
    if (isVector) {
        for (let i = 0; i < nPols; i++) {
            evals += `f${i+1}(𝔷),t${i+1}(𝔷),`;
        }
    }
    if (isSelected) evals += "fsel(𝔷),tsel(𝔷),";
    logger.info(`> STEP ${step}. Validate ${evals},S(𝔷·𝛚) ∈ 𝔽`);
    if(!validateEvaluations()) return false;
    ++step;

    let challs = "𝜸,𝜶,𝔷,v,u";
    if (isVector) challs = "𝛽," + challs;
    logger.info(`> STEP ${step}. Compute ${challs}`);
    computeChallenges();
    ++step;

    logger.info(`> STEP ${step}. Compute ZH(𝔷) and L₁(𝔷)`);
    const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
    const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);
    logger.info("··· ZH(𝔷) =", Fr.toString(ZHxi));
    logger.info("··· L₁(𝔷) =", Fr.toString(L1xi));
    ++step;

    logger.info(`> STEP ${step}. Compute r₀ = `);
    let r0 = Fr.zero;
    if (isSelected) {
        const selTBin = Fr.sub(proof.evaluations["selTxi"], Fr.square(proof.evaluations["selTxi"]));
        r0 = Fr.mul(Fr.add(r0, selTBin), challenges.alpha);

        const selFBin = Fr.sub(proof.evaluations["selFxi"], Fr.square(proof.evaluations["selFxi"]));
        r0 = Fr.mul(Fr.add(r0, selFBin), challenges.alpha);
    }

    let fxi = Fr.zero;
    let txi = Fr.zero;
    for (let i = nPols - 1; i >= 0; i--) {
        const nameEvalPolF = isVector ? `f${i}xi` : "fxi";
        const nameEvalPolT = isVector ? `t${i}xi` : "txi";

        fxi = Fr.add(Fr.mul(fxi, challenges.beta), proof.evaluations[nameEvalPolF]);
        txi = Fr.add(Fr.mul(txi, challenges.beta), proof.evaluations[nameEvalPolT]);
    }
    const fxigamma = Fr.add(fxi, challenges.gamma);
    const txigamma = Fr.add(txi, challenges.gamma);

    let r01 = Fr.mul(proof.evaluations["sxiw"], Fr.mul(fxigamma, txigamma));
    if (isSelected) {
        const selFGamma = Fr.mul(proof.evaluations["selFxi"], txigamma);
        const selTGamma = Fr.mul(proof.evaluations["selTxi"], fxigamma);

        r01 = Fr.add(r01, selTGamma);
        r01 = Fr.sub(r01, selFGamma);
    } else {
        r01 = Fr.add(r01, fxi);
        r01 = Fr.sub(r01, txi);
    }

    r0 = Fr.mul(Fr.add(r0, r01), challenges.alpha);
    logger.info("··· r₀    =", Fr.toString(r0));
    ++step;

    logger.info(`> STEP ${step}. Compute [D]₁ = `);
    let D1_12 = Fr.mul(challenges.alpha, fxigamma);
    D1_12 = Fr.mul(D1_12, txigamma);
    let D1_1 = Fr.add(Fr.sub(L1xi, D1_12), challenges.u);
    D1_1 = G1.timesFr(proof.commitments["S"], D1_1);
    const D1_2 = G1.timesFr(proof.commitments["Q"], ZHxi);
    const D1 = G1.sub(D1_1, D1_2);
    logger.info("··· [D]₁  =", G1.toString(G1.toAffine(D1)));
    ++step;

    logger.info(`> STEP ${step}. Compute [F]₁ = `);
    let F1 = G1.zero;
    if (isSelected) {
        F1 = G1.add(F1, proof.commitments["selT"]);
        F1 = G1.add(G1.timesFr(F1, challenges.v), proof.commitments["selF"]);
    }

    for (let i = nPols - 1; i >= 0; i--) {
        const namePolT = isVector ? `T${i}` : "T";

        F1 = G1.add(G1.timesFr(F1, challenges.v), proof.commitments[namePolT]);
    }
    for (let i = nPols - 1; i >= 0; i--) {
        const namePolF = isVector ? `F${i}` : "F";

        F1 = G1.add(G1.timesFr(F1, challenges.v), proof.commitments[namePolF]);
    }
    F1 = G1.add(G1.timesFr(F1, challenges.v), D1);
    logger.info("··· [F]₁  =", G1.toString(G1.toAffine(F1)));
    ++step;

    logger.info(`> STEP ${step}. Compute [E]₁ = `);
    let E1 = Fr.zero;
    if (isSelected) {
        E1 = Fr.add(E1, proof.evaluations["selTxi"]);
        E1 = Fr.add(Fr.mul(E1, challenges.v), proof.evaluations["selFxi"]);
    }

    for (let i = nPols - 1; i >= 0; i--) {
        const nameEvalPolT = isVector ? `t${i}xi` : "txi";

        E1 = Fr.add(Fr.mul(E1, challenges.v), proof.evaluations[nameEvalPolT]);
    }
    for (let i = nPols - 1; i >= 0; i--) {
        const nameEvalPolF = isVector ? `f${i}xi` : "fxi";

        E1 = Fr.add(Fr.mul(E1, challenges.v), proof.evaluations[nameEvalPolF]);
    }

    const E1_2 = Fr.mul(challenges.u, proof.evaluations["sxiw"]);
    E1 = Fr.add(Fr.mul(E1, challenges.v), E1_2);
    E1 = Fr.sub(E1, r0);
    E1 = G1.timesFr(G1.one, E1);
    logger.info("··· [E]₁  =", G1.toString(G1.toAffine(E1)));
    ++step;

    logger.info(`> STEP ${step}. Check pairing equation e(-[W𝔷(x)]₁ - u·[W𝔷·𝛚(x)]₁, [x]₂)·e(𝔷·[W𝔷(x)]₁ + u𝔷ω·[W𝔷·𝛚(x)]₁ + [F]₁ - [E]₁, [1]₂) = 1`);
    let A = G1.timesFr(proof.commitments["Wxiw"], challenges.u);
    A = G1.add(proof.commitments["Wxi"], A);

    let B = Fr.mul(challenges.u, Fr.w[nBits]);
    B = G1.timesFr(proof.commitments["Wxiw"], B);
    B = G1.add(proof.commitments["Wxi"], B);
    B = G1.timesFr(B, challenges.xi);
    B = G1.add(B, F1);
    B = G1.sub(B, E1);

    const isValid = await curve.pairingEq(G1.neg(A), X2, B, G2.one);

    if (isValid) {
        logger.info("> VERIFICATION OK");
    } else {
        logger.error("> VERIFICATION FAILED");
    }

    logger.info("> MULTISET EQUALITY KZG GRAND-SUM VERIFIER FINISHED");

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
        for (let i = 0; i < nPols; i++) {
            const namePolF = isVector ? `F${i}` : "F";
            const namePolT = isVector ? `T${i}` : "T";
            const lognamePolF = isVector ? `f${i+1}(x)` : "f(x)";
            const lognamePolT = isVector ? `t${i+1}(x)` : "t(x)";

            if (!valueBelongsToGroup1(`[${lognamePolF}]₁`, proof.commitments[namePolF])) return false;
            if (!valueBelongsToGroup1(`[${lognamePolT}]₁`, proof.commitments[namePolT])) return false;
        }

        if (isSelected) {
            if (!valueBelongsToGroup1("[fsel(x)]₁", proof.commitments["selF"])) return false;
            if (!valueBelongsToGroup1("[tsel(x)]₁", proof.commitments["selT"])) return false;
        }

        return (
            valueBelongsToGroup1("[S(x)]₁", proof.commitments["S"]) &&
            valueBelongsToGroup1("[Q(x)]₁", proof.commitments["Q"]) &&
            valueBelongsToGroup1("[W𝔷(x)]₁", proof.commitments["Wxi"]) &&
            valueBelongsToGroup1("[W𝔷·𝛚(x)]₁", proof.commitments["Wxiw"])
        );
    }

    function validateEvaluations() {
        for (let i = 0; i < nPols; i++) {
            const nameEvalPolF = isVector ? `f${i}xi` : "fxi";
            const nameEvalPolT = isVector ? `t${i}xi` : "txi";
            const lognameEvalPolF = isVector ? `f${i+1}(𝔷)` : "f(𝔷)";
            const lognameEvalPolT = isVector ? `t${i+1}(𝔷)` : "t(𝔷)";

            if (!valueBelongsToField(`${lognameEvalPolF}`, proof.evaluations[nameEvalPolF])) return false;
            if (!valueBelongsToField(`${lognameEvalPolT}`, proof.evaluations[nameEvalPolT])) return false;
        }

        return valueBelongsToField("S(𝔷·𝛚)", proof.evaluations["sxiw"]);
    }

    function computeChallenges() {
        const transcript = new Keccak256Transcript(curve);

        // STEP 1.1 Calculate challenge beta from transcript
        for (let i = 0; i < nPols; i++) {
            const namePolF = isVector ? `F${i}` : "F";
            const namePolT = isVector ? `T${i}` : "T";

            transcript.addPolCommitment(proof.commitments[namePolF]);
            transcript.addPolCommitment(proof.commitments[namePolT]);
        }

        if (isSelected) {
            transcript.addPolCommitment(proof.commitments["selF"]);
            transcript.addPolCommitment(proof.commitments["selT"]);
        }

        if (isVector) {
            challenges.beta = transcript.getChallenge();
            logger.info("··· 𝛃 =", Fr.toString(challenges.beta));

            transcript.addFieldElement(challenges.beta);
        }

        // STEP 1.2 Calculate challenge gamma from transcript
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
        for (let i = 0; i < nPols; i++) { // TODO (Héctor): Is this order the most appropriate one?
            const nameEvalPolF = isVector ? `f${i}xi` : "fxi";
            const nameEvalPolT = isVector ? `t${i}xi` : "txi";

            transcript.addFieldElement(proof.evaluations[nameEvalPolF]);
            transcript.addFieldElement(proof.evaluations[nameEvalPolT]);
        }

        if (isSelected) {
            transcript.addFieldElement(proof.evaluations["selFxi"]);
            transcript.addFieldElement(proof.evaluations["selTxi"]);
        }

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
};