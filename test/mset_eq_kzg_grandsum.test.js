const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const {
    getRandomValue,
    getBufferFromArray,
    getRandomBuffer,
    getRandomBinBuffer,
    getOneBuffer
} = require("./test.utils.js");
const path = require("path");

const mset_eq_kzg_grandsum_prover = require("../src/mset_eq_kzg_grandsum_prover.js");
const mset_eq_kzg_grandsum_verifier = require("../src/mset_eq_kzg_grandsum_verifier.js");
const { Polynomial } = require("../src/polynomial/polynomial.js");

describe("Protocols based on grand-sums", function () {
    this.timeout(1000000);

    let curve;

    before(async () => {
        curve = await getCurveFromName("bn128");
    });

    after(async () => {
        await curve.terminate();
    });

    it("Should perform the full proving and verifying process of a SNARK for multiset equalities, based on grand-sums and KZG", async () => {
        const nBits =  getRandomValue(1, 10);

        let evalsBufferF = new Array(1);
        let evalsBufferT = new Array(1);
        evalsBufferF[0] = getRandomBuffer(2 ** nBits, curve);
        evalsBufferT[0] = new Uint8Array(evalsBufferF[0].byteLength);

        evalsBufferT[0].set(evalsBufferF[0].slice(0, evalsBufferF[0].byteLength - curve.Fr.n8), curve.Fr.n8);
        evalsBufferT[0].set(evalsBufferF[0].slice(evalsBufferF[0].byteLength - curve.Fr.n8, evalsBufferF[0].byteLength), 0);

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT);

        const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        assert.ok(isValid);
    });

    it("Should perform the full proving and verifying process of a SNARK for vector multiset equalities, based on grand-sums and KZG", async () => {
        const nBits =  getRandomValue(1, 10);
        const nPols =  getRandomValue(2, 10);

        let evalsBufferF = new Array(nPols);
        let evalsBufferT = new Array(nPols);
        for (let i = 0; i < nPols; i++) {
            evalsBufferF[i] = getRandomBuffer(2 ** nBits, curve);
            evalsBufferT[i] = new Uint8Array(evalsBufferF[i].byteLength);

            evalsBufferT[i].set(evalsBufferF[i].slice(0, evalsBufferF[i].byteLength - curve.Fr.n8), curve.Fr.n8);
            evalsBufferT[i].set(evalsBufferF[i].slice(evalsBufferF[i].byteLength - curve.Fr.n8, evalsBufferF[i].byteLength), 0);
        }

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT);

        const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        assert.ok(isValid);
    });

    it("Should perform the full proving and verifying process of a SNARK for a simple selected vector multiset equalities, based on grand-sums and KZG", async () => {
        const nBits =  getRandomValue(2, 2);
        const nPols =  1

        let one = curve.Fr.one;
        let two = curve.Fr.add(curve.Fr.one, curve.Fr.one);
        let three = curve.Fr.add(curve.Fr.add(curve.Fr.one, curve.Fr.one), curve.Fr.one);
        let four = curve.Fr.add(curve.Fr.add(curve.Fr.add(curve.Fr.one, curve.Fr.one), curve.Fr.one), curve.Fr.one);

        let testArray = [one, two, three, four];
        let evalsBufferF = new Array(nPols);
        let evalsBufferT = new Array(nPols);
        for (let i = 0; i < nPols; i++) {
            evalsBufferF[i] = getBufferFromArray(testArray, curve);
            evalsBufferT[i] = new Uint8Array(evalsBufferF[i].byteLength);

            evalsBufferT[i].set(evalsBufferF[i].slice(0, evalsBufferF[i].byteLength - curve.Fr.n8), curve.Fr.n8);
            evalsBufferT[i].set(evalsBufferF[i].slice(evalsBufferF[i].byteLength - curve.Fr.n8, evalsBufferF[i].byteLength), 0);
        }

        for (let i = 0; i < 4; i++) {
            console.log("evalsBufferF[i]: ", curve.Fr.toString(evalsBufferF[0].slice(i * curve.Fr.n8, i * curve.Fr.n8 + curve.Fr.n8)));
            console.log("evalsBufferT[i]: ", curve.Fr.toString(evalsBufferT[0].slice(i * curve.Fr.n8, i * curve.Fr.n8 + curve.Fr.n8)));
        }


        // Make it simple for now
        let evalsBufferFSelected = getOneBuffer(2 ** nBits, curve)
        let evalsBufferTSelected = getOneBuffer(2 ** nBits, curve)
        evalsBufferFSelected.set(curve.Fr.zero, 0);
        evalsBufferTSelected.set(curve.Fr.zero, curve.Fr.n8);

        for (let i = 0; i < 4; i++) {
            console.log("evalsBufferFSelected[i]: ", curve.Fr.toString(evalsBufferFSelected.slice(i * curve.Fr.n8, i * curve.Fr.n8 + curve.Fr.n8)));
            console.log("evalsBufferTSelected[i]: ", curve.Fr.toString(evalsBufferTSelected.slice(i * curve.Fr.n8, i * curve.Fr.n8 + curve.Fr.n8)));
        }

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT, evalsBufferFSelected, evalsBufferTSelected);

        const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
        assert.ok(isValid);
    });

    // it("Should perform the full proving and verifying process of a SNARK for selected vector multiset equalities, based on grand-sums and KZG", async () => {
    //     const nBits =  getRandomValue(2, 2);
    //     const nPols =  1

    //     let evalsBufferF = new Array(nPols);
    //     let evalsBufferT = new Array(nPols);
    //     for (let i = 0; i < nPols; i++) {
    //         evalsBufferF[i] = getRandomBuffer(2 ** nBits, curve);
    //         evalsBufferT[i] = new Uint8Array(evalsBufferF[i].byteLength);

    //         evalsBufferT[i].set(evalsBufferF[i].slice(0, evalsBufferF[i].byteLength - curve.Fr.n8), curve.Fr.n8);
    //         evalsBufferT[i].set(evalsBufferF[i].slice(evalsBufferF[i].byteLength - curve.Fr.n8, evalsBufferF[i].byteLength), 0);
    //     }

    //     // Make it simple for now
    //     let evalsBufferFSelected = getOneBuffer(2 ** nBits, curve)
    //     let evalsBufferTSelected = getOneBuffer(2 ** nBits, curve)
    //     evalsBufferFSelected.set(curve.Fr.zero, 0);
    //     evalsBufferTSelected.set(curve.Fr.zero, curve.Fr.n8);

    //     const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
    //     const proof = await mset_eq_kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT, evalsBufferFSelected, evalsBufferTSelected);

    //     const isValid = await mset_eq_kzg_grandsum_verifier(pTauFilename, proof, nBits);
    //     assert.ok(isValid);
    // });
});