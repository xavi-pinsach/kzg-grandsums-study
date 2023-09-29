const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("../Keccak256Transcript");
const { Polynomial } = require("../polynomial/polynomial");
const { Evaluations } = require("../polynomial/evaluations");
const { computeZHEvaluation, computeL1Evaluation } = require("../polynomial/polynomial_utils");
const readPTauHeader = require("../ptau_utils");
const ComputeZGrandProductPolynomial = require("./grandproduct");

const logger = require("../../logger.js");

module.exports = async function mset_eq_kzg_grandproduct_prover(pTauFilename, evalsFs, evalsTs) {
    logger.info("> MULTISET EQUALITY KZG GRAND-PRODUCT PROVER STARTED");

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve, power: nBitsPTau } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const sG1 = G1.F.n8 * 2;

    // The following are done to avoid the user having to provide input buffers of one polynomial as an array one element
    if (!Array.isArray(evalsFs)) {
        evalsFs = [evalsFs];
    }
    if (!Array.isArray(evalsTs)) {
        evalsTs = [evalsTs];
    }

    // Sanity checks
    if (evalsFs.length !== evalsTs.length) {
        throw new Error(`The lengths of the two vector multisets must be the same.`);
    }
    const nPols = evalsFs.length;
    if (nPols === 0) {
        throw new Error(`The number of multisets must be greater than 0.`);
    }

    // Ensure all polynomials have the same length
    for (let i = 0; i < nPols; i++) {
        if (evalsFs[i].length() !== evalsTs[i].length()) {
            throw new Error(`The ${i}-th multiset buffers must have the same length.`);
        } else if (evalsFs[i].length() !== evalsFs[0].length()) {
            throw new Error("The multiset buffers must all have the same length.");
        }
    }

    const nBits = Math.ceil(Math.log2(evalsFs[0].length()));
    const domainSize = 2 ** nBits;

    // Ensure the polynomial has a length that is equal to a power of two
    if (evalsFs[0].length() !== domainSize) {
        throw new Error("Polynomial length must be a power of two.");
    }

    // Ensure the powers of Tau file is sufficiently large
    if (nBitsPTau < nBits) {
        throw new Error("The Powers of Tau file is not sufficiently large to commit the polynomials.");
    }

    const PTau = new BigBuffer(domainSize * 2 * sG1);
    await fdPTau.readToBuffer(PTau, 0, domainSize * 2 * sG1, pTauSections[2][0].p);
    await fdPTau.close();

    logger.info("-------------------------------------");
    logger.info("  MULTISET EQUALITY KZG GRAND-PRODUCT PROVER SETTINGS");
    logger.info(`  Curve:       ${curve.name}`);
    logger.info(`  Domain size: ${domainSize}`);
    logger.info(`  Number of polynomials: ${nPols}`);
    // logger.info(`  Selectors: ${isSelected ? "Yes" : "No"}`);
    logger.info("-------------------------------------");

    let proof = {evaluations: {}, commitments: {}};
    let challenges = {};
    let polFs = new Array(nPols);
    let polTs = new Array(nPols);
    let polF, polT, evalsF, evalsT;
    let selF, selT;
    let polS;
    let polQ = Polynomial.zero(domainSize, curve);
    let polWxi = Polynomial.zero(domainSize, curve) 
    let polWxiomega = Polynomial.zero(domainSize, curve)

    const transcript = new Keccak256Transcript(curve);

    const isVector = nPols > 1;
    let round = 1;

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

    logger.info("> MULTISET EQUALITY KZG GRAND-PRODUCT PROVER FINISHED");

    return proof;

    async function computeWitnessPolynomials() {
        for (let i = 0; i < nPols; i++) {
            // Convert the evaluations to Montgomery form
            evalsFs[i].eval = await Fr.batchToMontgomery(evalsFs[i].eval);
            evalsTs[i].eval = await Fr.batchToMontgomery(evalsTs[i].eval);

            // Get the polynomials from the evaluations
            polFs[i] = await Polynomial.fromEvaluations(evalsFs[i].eval, curve);
            polTs[i] = await Polynomial.fromEvaluations(evalsTs[i].eval, curve);
        }

        for (let i = 0; i < nPols; i++) {
            const namePolF = isVector ? `F${i}` : "F";
            const namePolT = isVector ? `T${i}` : "T";
            const lognamePolF = isVector ? `f${i+1}(x)` : "f(x)";
            const lognamePolT = isVector ? `t${i+1}(x)` : "t(x)";

            proof.commitments[namePolF] = await commit(polFs[i]);
            proof.commitments[namePolT] = await commit(polTs[i]);

            logger.info(`··· [${lognamePolF}]₁ =`, G1.toString(proof.commitments[namePolF]));
            logger.info(`··· [${lognamePolT}]₁ =`, G1.toString(proof.commitments[namePolT]));
        }
    }

    async function ComputeZPolynomial() {
        // 1. Compute the challenges
        for (let i = 0; i < nPols; i++) {
            const namePolF = isVector ? `F${i}` : "F";
            const namePolT = isVector ? `T${i}` : "T";

            transcript.addPolCommitment(proof.commitments[namePolF]);
            transcript.addPolCommitment(proof.commitments[namePolT]);
        }

        if (isVector) {
            challenges.beta = transcript.getChallenge();
            logger.info("···      𝛃  =", Fr.toString(challenges.beta));

            transcript.addFieldElement(challenges.beta);
        }

        challenges.gamma = transcript.getChallenge();
        logger.info("···      𝜸  =", Fr.toString(challenges.gamma));

        // 2. Compute the random linear combination of the polynomials fᵢ,tᵢ ∈ 𝔽[X]
        if (isVector) {
            // QUESTION (Héctor): Should I compute them directly from the evaluations?
            polF = Polynomial.zero(domainSize, curve);
            polT = Polynomial.zero(domainSize, curve);
            for (let i = nPols - 1; i >= 0; i--) {
                polF.mulScalar(challenges.beta).add(polFs[i]);
                polT.mulScalar(challenges.beta).add(polTs[i]);
            }

            evalsF = await Evaluations.fromPolynomial(polF, 1, curve);
            evalsT = await Evaluations.fromPolynomial(polT, 1, curve);
        } else {
            polF = polFs[0];
            polT = polTs[0];

            evalsF = evalsFs[0];
            evalsT = evalsTs[0];
        }
    
        polZ = await ComputeZGrandProductPolynomial(evalsF, evalsT, challenges.gamma, curve);
    
        proof.commitments["Z"] = await commit(polZ);
        logger.info(`··· [Z(x)]₁ =`, G1.toString(proof.commitments["Z"]));
    }

    async function computeQPolynomial() {
        transcript.addFieldElement(challenges.gamma);
        transcript.addPolCommitment(proof.commitments["Z"]);

        challenges.alpha = transcript.getChallenge();
        logger.info("···      𝜶  =", Fr.toString(challenges.alpha));

        const polZ1 = polZ.clone();
        const polL1 = await Polynomial.Lagrange1(nBits, curve);
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