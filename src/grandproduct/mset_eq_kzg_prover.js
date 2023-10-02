const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("../Keccak256Transcript");
const { Polynomial } = require("../polynomial/polynomial");
const { Evaluations } = require("../polynomial/evaluations");
const { computeZHEvaluation, computeL1Evaluation } = require("../polynomial/polynomial_utils");
const { readPTauHeader } = require("../ptau_utils");
const { ComputeZGrandProductPolynomial } = require("./grandproduct");

const logger = require("../../logger.js");

module.exports = async function mset_eq_kzg_grandproduct_prover(pTauFilename, evalsFs, evalsTs, evalsSelF = null, evalsSelT = null) {
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

    // If the selection buffers are not provided, assume all elements are selected
    if (evalsSelF === null) {
        evalsSelF = Evaluations.getOneEvals(evalsFs[0].length(), curve);
    } 
    if (evalsSelT === null) {
        evalsSelT = Evaluations.getOneEvals(evalsTs[0].length(), curve);
    }

    // More sanity checks
    if (evalsSelF.length() !== evalsSelT.length()) {
        throw new Error("The selection buffers must have the same length.");
    } else if (evalsSelF.length() !== evalsFs[0].length()) {
        throw new Error("The selection buffers must have the same length as the multiset buffers.");
    }

    // Checking for "trivial" cases
    let isSelected = true;
    if (evalsSelF.isAllOnes() && evalsSelT.isAllOnes()) {
        isSelected = false;
    } else if (evalsSelF.isAllZeros() && evalsSelT.isAllZeros()) {
        logger.warn("The selection buffers are all zeros. The argument is trivially satisfied.");
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
    logger.info(`  Selectors: ${isSelected ? "Yes" : "No"}`);
    logger.info("-------------------------------------");

    let proof = {evaluations: {}, commitments: {}};
    let challenges = {};
    let polFs = new Array(nPols);
    let polTs = new Array(nPols);
    let polF, polT, evalsF, evalsT;
    let selF, selT;
    let polZ;
    let polQ = Polynomial.zero(domainSize, curve);
    let polWxi = Polynomial.zero(domainSize, curve) 
    let polWxiomega = Polynomial.zero(domainSize, curve)

    const transcript = new Keccak256Transcript(curve);

    const isVector = nPols > 1;
    let round = 1;

    let msg = "> ROUND ${round}. Generate the witness polynomials";
    if (isVector) {
        msg += ` fᵢ,tᵢ ∈ 𝔽[X], for i ∈ [${nPols}]`;
    } else {
        msg += ` f,t ∈ 𝔽[X]`;
    }
    if (isSelected) {
        msg += ", and the selector polynomials fsel,tsel ∈ 𝔽[X]";
    }
    logger.info(msg);
    await computeWitnessPolynomials();
    ++round;

    logger.info(`> ROUND ${round}. Compute the grand-product polynomial Z ∈ 𝔽[X]`);
    await ComputeZPolynomial();
    ++round;

    logger.info(`> ROUND ${round}. Compute the quotient polynomial Q ∈ 𝔽[X]`);
    await computeQPolynomial();
    ++round;

    logger.info(`> ROUND ${round}. Compute the evaluations of the polynomials`);
    computeEvaluations();
    ++round;

    logger.info(`> ROUND ${round}. Compute the opening proof polynomials W𝔷, W𝔷𝛚 ∈ 𝔽[X]`);
    await computeW();

    logger.info("");
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

        // Compute the selector polynomials
        if (isSelected) {
            selF = await Polynomial.fromEvaluations(evalsSelF.eval, curve);
            selT = await Polynomial.fromEvaluations(evalsSelT.eval, curve);

            proof.commitments["selF"] = await commit(selF);
            proof.commitments["selT"] = await commit(selT);
            
            logger.info(`··· [fsel(x)]₁ =`, G1.toString(proof.commitments["selF"]));
            logger.info(`··· [tsel(x)]₁ =`, G1.toString(proof.commitments["selT"]));
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

        if (isSelected) {
            transcript.addPolCommitment(proof.commitments["selF"]);
            transcript.addPolCommitment(proof.commitments["selT"]);
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
    
        // 3. Compute the grand-product polynomial S(X)
        polZ = await ComputeZGrandProductPolynomial(evalsF, evalsT, evalsSelF, evalsSelT, isSelected, challenges.gamma, curve);

        proof.commitments["Z"] = await commit(polZ);
        logger.info(`··· [Z(x)]₁ =`, G1.toString(proof.commitments["Z"]));
    }

    async function computeQPolynomial() {
        transcript.addFieldElement(challenges.gamma);
        transcript.addPolCommitment(proof.commitments["Z"]);

        challenges.alpha = transcript.getChallenge();
        logger.info("···      𝜶  =", Fr.toString(challenges.alpha));

        if (isSelected) {
            const selTBin1 = selT.clone()
            await selTBin1.multiply(selT.clone());
            const selTBin = selT.clone().sub(selTBin1);
            polQ.add(selTBin).mulScalar(challenges.alpha);

            const selFBin1 = selF.clone()
            await selFBin1.multiply(selF.clone());
            const selFBin = selF.clone().sub(selFBin1);
            polQ.add(selFBin).mulScalar(challenges.alpha);
        }

        const polQ1 = polZ.clone();
        await polQ1.shiftOmega();
        const polQ2 = polZ.clone();

        const polFGamma = polF.clone().addScalar(challenges.gamma);
        const polTGamma = polT.clone().addScalar(challenges.gamma);

        if (isSelected) {
            polTGamma.subScalar(Fr.one);
            await polTGamma.multiply(selT.clone());
            polTGamma.addScalar(Fr.one);
            await polQ1.multiply(polTGamma);

            polFGamma.subScalar(Fr.one);
            await polFGamma.multiply(selF.clone());
            polFGamma.addScalar(Fr.one);
            await polQ2.multiply(polFGamma);
        } else {
            await polQ1.multiply(polTGamma);
            await polQ2.multiply(polFGamma);
        }

        polQ1.sub(polQ2);
        polQ.add(polQ1).mulScalar(challenges.alpha);

        const polQ3 = polZ.clone();
        const polL1 = await Polynomial.Lagrange1(nBits, curve);
        polQ3.subScalar(Fr.one);
        await polQ3.multiply(polL1);
        polQ.add(polQ3);

        polQ.divZh(domainSize);

        proof.commitments["Q"] = await commit(polQ);
        logger.info(`··· [Q(x)]₁ =`, G1.toString(proof.commitments["Q"]));
    }

    function computeEvaluations() {
        transcript.addFieldElement(challenges.alpha);
        transcript.addPolCommitment(proof.commitments["Q"]);

        challenges.xi = transcript.getChallenge();
        logger.info("···      𝔷  =", Fr.toString(challenges.xi));

        for (let i = 0; i < nPols; i++) {
            const nameEvalPolF = isVector ? `f${i}xi` : "fxi";
            const lognameEvalPolF = isVector ? `f${i+1}(𝔷)` : "f(𝔷)";

            proof.evaluations[nameEvalPolF] = polFs[i].evaluate(challenges.xi);

            logger.info(`···   ${lognameEvalPolF}  =`, Fr.toString(proof.evaluations[nameEvalPolF]));
        }

        if (isSelected) {
            proof.evaluations["selFxi"] = selF.evaluate(challenges.xi);
            proof.evaluations["selTxi"] = selT.evaluate(challenges.xi);

            logger.info(`···   fsel(𝔷)  =`, Fr.toString(proof.evaluations["selFxi"]));
            logger.info(`···   tsel(𝔷)  =`, Fr.toString(proof.evaluations["selTxi"]));
        }

        proof.evaluations["zxiw"] = polZ.evaluate(Fr.mul(challenges.xi, Fr.w[nBits]));
        logger.info(`··· Z(𝔷·𝛚)  =`, Fr.toString(proof.evaluations["zxiw"]));
    }

    async function computeW() {
        transcript.addFieldElement(challenges.xi);
        for (let i = 0; i < nPols; i++) {
            const nameEvalPolF = isVector ? `f${i}xi` : "fxi";
            transcript.addFieldElement(proof.evaluations[nameEvalPolF]);
        }

        if (isSelected) {
            transcript.addFieldElement(proof.evaluations["selFxi"]);
            transcript.addFieldElement(proof.evaluations["selTxi"]);
        }

        transcript.addFieldElement(proof.evaluations["zxiw"]);

        challenges.v = transcript.getChallenge();
        logger.info("···      v  = ", Fr.toString(challenges.v));

        const ZHxi = computeZHEvaluation(curve, challenges.xi, nBits);
        const L1xi = computeL1Evaluation(curve, challenges.xi, ZHxi, nBits);

        logger.info("···  ZH(𝔷)  =", Fr.toString(ZHxi));
        logger.info("···  L₁(𝔷)  =", Fr.toString(L1xi));

        // Compute the linearisation polynomial r(X)
        let polR = Polynomial.zero(domainSize, curve);
        if (isSelected) {
            const selTBin = Fr.sub(proof.evaluations["selTxi"], Fr.square(proof.evaluations["selTxi"]));
            polR.addScalar(selTBin).mulScalar(challenges.alpha);

            const selFBin = Fr.sub(proof.evaluations["selFxi"], Fr.square(proof.evaluations["selFxi"]));
            polR.addScalar(selFBin).mulScalar(challenges.alpha);
        }

        const polR1 = Polynomial.zero(domainSize, curve);

        const fxi = polF.evaluate(challenges.xi);
        let fxigamma = Fr.add(fxi, challenges.gamma);
        let txigamma = polT.addScalar(challenges.gamma);

        if (isSelected) {
            fxigamma = Fr.sub(fxigamma, Fr.one);
            txigamma = txigamma.subScalar(Fr.one);
            const selFGamma = Fr.add(Fr.mul(proof.evaluations["selFxi"], fxigamma), Fr.one);
            const selTGamma = txigamma.mulScalar(proof.evaluations["selTxi"]).addScalar(Fr.one);

            selTGamma.mulScalar(proof.evaluations["zxiw"]);
            polR1.add(selTGamma);

            const ZZ = polZ.clone();
            ZZ.mulScalar(selFGamma);

            polR1.sub(ZZ);
        } else {
            txigamma.mulScalar(proof.evaluations["zxiw"]);
            polR1.add(txigamma);

            const ZZ = polZ.clone();
            ZZ.mulScalar(fxigamma);

            polR1.sub(ZZ);
        }

        polR.add(polR1).mulScalar(challenges.alpha);

        const polR2 = polZ.clone().subScalar(Fr.one);
        polR2.mulScalar(L1xi);
        polR.add(polR2);

        const polR3 = polQ.clone().mulScalar(ZHxi);
        polR.sub(polR3);

        // Compute the polynomial W𝔷(X)
        if (isSelected) {
            polWxi.add(selT.clone().subScalar(proof.evaluations["selTxi"]));
            polWxi.mulScalar(challenges.v).add(selF.clone().subScalar(proof.evaluations["selFxi"]));
        }

        for (let i = nPols - 1; i >= 0; i--) {
            const nameEvalPolF = isVector ? `f${i}xi` : "fxi";
            polWxi.mulScalar(challenges.v).add(polFs[i].clone().subScalar(proof.evaluations[nameEvalPolF]));
        }

        polWxi.mulScalar(challenges.v).add(polR.clone());
        polWxi.divByXSubValue(challenges.xi);

        // Compute the polynomial W𝔷·𝛚(X)
        polWxiomega = polZ.clone().subScalar(proof.evaluations["zxiw"]);
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