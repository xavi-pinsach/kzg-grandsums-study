const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const { getRandomValue } = require("./test.utils.js");
const path = require("path");

const mset_eq_kzg_grandsum_prover = require("../src/mset_eq_kzg_grandsum_prover.js");
const mset_eq_kzg_grandsum_verifier = require("../src/mset_eq_kzg_grandsum_verifier.js");
const { Evaluations } = require("../src/polynomial/evaluations.js");

describe("Protocols based on grand-sums and KZG", () => {
    let curve;
    let Fr;

    before(async () => {
        curve = await getCurveFromName("bn128");
        Fr = curve.Fr;
    });

    after(async () => {
        await curve.terminate();
    });

    describe("Full proving and verifying process of a SNARK for multiset equalities", () => {
        it("Should proof and verify a standard multiset equality", async () => {
            const nBits =  getRandomValue(1, 10);

            let evalsF = Evaluations.getRandomEvals(2 ** nBits, curve);
            let evalsT = Evaluations.fromEvals(evalsF);
            evalsT.setEvaluation(1, evalsF.getEvaluationSequence(0, evalsF.length() - 1));
            evalsT.setEvaluation(0, evalsF.getEvaluation(evalsF.length() - 1));

            const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
            const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT);

            const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
            assert.ok(isValid);
        });

        it("Should proof and verify a vector multiset equality", async () => {
            const nBits =  getRandomValue(1, 10);
            const nPols =  getRandomValue(2, 10);

            let evalsF = new Array(nPols);
            let evalsT = new Array(nPols);
            for (let i = 0; i < nPols; i++) {
                evalsF[i] = Evaluations.getRandomEvals(2 ** nBits, curve);
                evalsT[i] = Evaluations.fromEvals(evalsF[i])
                evalsT[i].setEvaluation(1, evalsF[i].getEvaluationSequence(0, evalsF[i].length() - 1));
                evalsT[i].setEvaluation(0, evalsF[i].getEvaluation(evalsF[i].length() - 1));
            }

            const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
            const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT);

            const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
            assert.ok(isValid);
        });

        it("Should proof and verify a selected multiset equality", async () => {
            const nBits =  getRandomValue(1, 10);

            evalsF = Evaluations.getRandomEvals(2 ** nBits, curve);
            evalsT = Evaluations.fromEvals(evalsF)
            evalsT.setEvaluation(1, evalsF.getEvaluationSequence(0, evalsF.length() - 1));
            evalsT.setEvaluation(0, evalsF.getEvaluation(evalsF.length() - 1));

            // Make the selection simple (but non trivial) for now
            let evalsFSelected = Evaluations.getOneEvals(2 ** nBits, curve)
            let evalsTSelected = Evaluations.fromEvals(evalsFSelected)
            evalsFSelected.setEvaluation(evalsFSelected.length() - 1, Fr.zero);
            evalsTSelected.setEvaluation(0, Fr.zero);

            const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
            const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT, evalsFSelected, evalsTSelected);

            const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
            assert.ok(isValid);
        });

        it("Should proof and verify a selected vector multiset equality", async () => {
            const nBits =  getRandomValue(1, 10);
            const nPols =  getRandomValue(2, 10);

            let evalsF = new Array(nPols);
            let evalsT = new Array(nPols);
            for (let i = 0; i < nPols; i++) {
                evalsF[i] = Evaluations.getRandomEvals(2 ** nBits, curve);
                evalsT[i] = Evaluations.fromEvals(evalsF[i])
                evalsT[i].setEvaluation(1, evalsF[i].getEvaluationSequence(0, evalsF[i].length() - 1));
                evalsT[i].setEvaluation(0, evalsF[i].getEvaluation(evalsF[i].length() - 1));
            }

            // Make the selection simple (but non trivial) for now
            let evalsFSelected = Evaluations.getOneEvals(2 ** nBits, curve)
            let evalsTSelected = Evaluations.fromEvals(evalsFSelected)
            evalsFSelected.setEvaluation(evalsFSelected.length() - 1, Fr.zero);
            evalsTSelected.setEvaluation(0, Fr.zero);

            const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
            const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT, evalsFSelected, evalsTSelected);

            const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
            assert.ok(isValid);
        });
    });
});