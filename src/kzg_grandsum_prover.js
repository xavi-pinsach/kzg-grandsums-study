const { readBinFile } = require("@iden3/binfileutils");
const { BigBuffer } = require("ffjavascript");
const { Keccak256Transcript } = require("./Keccak256Transcript");
const { Polynomial } = require("./polynomial/polynomial");
const { Evaluations } = require("./polynomial/evaluations");
const ComputeSGrandSumPolynomial = require("./grandsum");
const readPTauHeader = require("./ptau_utils");
const { computeZHEvaluation, computeL1Evaluation } = require("./polynomial/polynomial_utils");

const logger = require("../logger.js");

module.exports = async function kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT, nPols = 1) {
    logger.info("> KZG GRAND SUM PROVER STARTED");

    if (nPols < 1) throw new Error("The number of polynomials must be greater than 0.");

    const { fd: fdPTau, sections: pTauSections } = await readBinFile(pTauFilename, "ptau", 1, 1 << 22, 1 << 24);
    const { curve, power: nBitsPTau } = await readPTauHeader(fdPTau, pTauSections);
    const Fr = curve.Fr;
    const G1 = curve.G1;
    const sG1 = G1.F.n8 * 2;

    let evalsFs = [];
    let evalsTs = [];
    for (let i = 0; i < nPols; i++) {
        evalsFs.push(new Evaluations(evalsBufferF[i], curve, logger));
        evalsTs.push(new Evaluations(evalsBufferT[i], curve, logger));

        // Ensure all polynomials have the same length
        if (evalsFs[i].length() !== evalsTs[i].length()) {
            throw new Error(`The ${i}-th buffers must have the same length.`);
        } else if (evalsFs[i].length() !== evalsFs[0].length()) {
            throw new Error("The buffers must all have the same length.");
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
    logger.info("  KZG GRAND SUM PROVER SETTINGS");
    logger.info(`  Curve:       ${curve.name}`);
    logger.info(`  Domain size: ${domainSize}`);
    logger.info("-------------------------------------");

    let proof = {evaluations: {}, commitments: {}};
    let challenges = {};
    let polFs = new Array(nPols);
    let polTs = new Array(nPols);
    let polF, polT, evalsF, evalsT;
    let polS, polQ;

    const transcript = new Keccak256Transcript(curve);

    const isVector = nPols > 1;
    let round = 1;
    if (isVector) {
        logger.info(`> ROUND ${round}. Generate the witness polynomials fᵢ,tᵢ ∈ 𝔽[X] from the evaluations, for i ∈ [${nPols}]`);
    } else {
        logger.info(`> ROUND ${round}. Generate the witness polynomials f,t ∈ 𝔽[X] from the evaluations`);
    }
    await computeWitnessPolynomials();
    ++round;

    if (isVector) {
        logger.info(`> ROUND ${round}. Generate the randomly combined polynomials f,t ∈ 𝔽[X]`);
        await computeRandCombPolynomials();
        ++round;
    }

    logger.info(`> ROUND ${round}. Compute the grand-sum polynomial S ∈ 𝔽[X]`);
    await ComputeSPolynomial();
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
    logger.info("> KZG GRAND SUM PROVER FINISHED");

    return proof;

    async function computeWitnessPolynomials() {
        for (let i = 0; i < nPols; i++) {
            // Convert the evaluations to Montgomery form
            evalsFs[i].eval = await Fr.batchToMontgomery(evalsFs[i].eval);
            evalsTs[i].eval = await Fr.batchToMontgomery(evalsTs[i].eval);

            // Get the polynomials from the evaluations
            polFs[i] = await Polynomial.fromEvaluations(evalsFs[i].eval, curve, logger);
            polTs[i] = await Polynomial.fromEvaluations(evalsTs[i].eval, curve, logger);

            // TODO: I am checking the vector condition at each iteration, can I do it only once without coding overload?
            if (isVector) {
                proof.commitments[`F${i}`] = await commit(polFs[i]);
                proof.commitments[`T${i}`] = await commit(polTs[i]);
    
                logger.info(`··· [f${i+1}(x)]₁ =`, G1.toString(proof.commitments[`F${i}`]));
                logger.info(`··· [t${i+1}(x)]₁ =`, G1.toString(proof.commitments[`T${i}`]));
            } else {
                evalsF = evalsFs[0];
                evalsT = evalsTs[0];
                polF = polFs[0];
                polT = polTs[0];

                proof.commitments["F"] = await commit(polF);
                proof.commitments["T"] = await commit(polT);
    
                logger.info("··· [f(x)]₁ =", G1.toString(proof.commitments["F"]));
                logger.info("··· [t(x)]₁ =", G1.toString(proof.commitments["T"]));
            }

        }
    }

    async function computeRandCombPolynomials() {
        for (let i = 0; i < nPols; i++) {
            transcript.addPolCommitment(proof.commitments[`F${i}`]);
            transcript.addPolCommitment(proof.commitments[`T${i}`]);
        }

        challenges.beta = transcript.getChallenge();
        logger.info("···      𝛃  =", Fr.toString(challenges.beta));

        polF = Polynomial.zero(curve, logger);
        polT = Polynomial.zero(curve, logger);
        for (let i = nPols - 1; i >= 0; i--) {
            polF.mulScalar(challenges.beta).add(polFs[i]);
            polT.mulScalar(challenges.beta).add(polTs[i]);
        }

        evalsF = await Evaluations.fromPolynomial(polF, 1, curve, logger);
        evalsT = await Evaluations.fromPolynomial(polT, 1, curve, logger);

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