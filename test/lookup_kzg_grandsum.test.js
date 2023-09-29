const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const { getRandomValue, getRandomBuffer } = require("./test.utils.js");
const { Evaluations } = require("../src/polynomial/evaluations");
const path = require("path");

const mset_eq_kzg_grandsum_prover = require("../src/grandsum/mset_eq_kzg_prover.js");
const mset_eq_kzg_grandsum_verifier = require("../src/grandsum/mset_eq_kzg_verifier.js");

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

    describe("Full proving and verifying process of a SNARK for lookups", () => {
        // it("Should proof and verify a standard lookup", async () => {
        //     const nBits =  getRandomValue(2,2);

        //     let evalsT = Evaluations.getRandomEvals(2 ** nBits, curve);
        //     let evalsF = Evaluations.fromEvals(evalsT);
        //     // Copy two elements from F to another position
        //     evalsF.setEvaluation(1, evalsF.getEvaluation(0));
        //     evalsF.setEvaluation(evalsF.length() - 1, evalsF.getEvaluation(0));

        //     let evalsMulCounter = Evaluations.getOneEvals(2 ** nBits, curve);
        //     evalsMulCounter.setEvaluation(0, Fr.e(3));
        //     evalsMulCounter.setEvaluation(1, Fr.zero);
        //     evalsMulCounter.setEvaluation(evalsMulCounter.length() - 1, Fr.zero);

        //     const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        //     // TODO: Refactor to accept Evaluations
        //     const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF.eval.slice(), evalsT.eval.slice(), Evaluations.getOneEvals(2 ** nBits, curve).eval.slice(), evalsMulCounter.eval.slice());

        //     const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        //     assert.ok(isValid);
        // });

        // it("Should proof and verify a vector lookup", async () => {
        //     const nBits =  getRandomValue(1, 10);
        //     const nPols =  getRandomValue(2, 10);

        //     let evalsF = new Array(nPols);
        //     let evalsT = new Array(nPols);
        //     for (let i = 0; i < nPols; i++) {
        //         evalsF[i] = getRandomBuffer(2 ** nBits, curve);
        //         evalsT[i] = new Uint8Array(evalsF[i].byteLength);
        //         evalsT[i].set(evalsF[i].slice(0, evalsF[i].byteLength - Fr.n8), Fr.n8);
        //         evalsT[i].set(evalsF[i].slice(evalsF[i].byteLength - Fr.n8, evalsF[i].byteLength), 0);
        //     }

        //     const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        //     const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT);

        //     const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        //     assert.ok(isValid);
        // });

        // it("Should proof and verify a selected lookup", async () => {
        //     const nBits =  getRandomValue(1, 10);

        //     evalsF = getRandomBuffer(2 ** nBits, curve);
        //     evalsT = new Uint8Array(evalsF.byteLength);
        //     evalsT.set(evalsF.slice(0, evalsF.byteLength - Fr.n8), Fr.n8);
        //     evalsT.set(evalsF.slice(evalsF.byteLength - Fr.n8, evalsF.byteLength), 0);

        //     // Make the selection simple (but non trivial) for now
        //     let evalsFSelected = getOneBuffer(2 ** nBits, curve)
        //     let evalsTSelected = getOneBuffer(2 ** nBits, curve)
        //     evalsFSelected.set(Fr.zero, 0);
        //     evalsTSelected.set(Fr.zero, Fr.n8);

        //     const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        //     const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT, evalsFSelected, evalsTSelected);

        //     const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        //     assert.ok(isValid);
        // });

        // it("Should proof and verify a selected vector lookup", async () => {
        //     const nBits =  getRandomValue(1, 10);
        //     const nPols =  getRandomValue(2, 10);

        //     let evalsF = new Array(nPols);
        //     let evalsT = new Array(nPols);
        //     for (let i = 0; i < nPols; i++) {
        //         evalsF[i] = getRandomBuffer(2 ** nBits, curve);
        //         evalsT[i] = new Uint8Array(evalsF[i].byteLength);
        //         evalsT[i].set(evalsF[i].slice(0, evalsF[i].byteLength - Fr.n8), Fr.n8);
        //         evalsT[i].set(evalsF[i].slice(evalsF[i].byteLength - Fr.n8, evalsF[i].byteLength), 0);
        //     }

        //     // Make the selection simple (but non trivial) for now
        //     let evalsFSelected = getOneBuffer(2 ** nBits, curve)
        //     let evalsTSelected = getOneBuffer(2 ** nBits, curve)
        //     evalsFSelected.set(Fr.zero, 0);
        //     evalsTSelected.set(Fr.zero, Fr.n8);

        //     const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        //     const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsF, evalsT, evalsFSelected, evalsTSelected);

        //     const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        //     assert.ok(isValid);
        // });
    });
});