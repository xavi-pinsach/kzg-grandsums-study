const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const {
    getRandomValue,
    getRandomBuffer,
} = require("./test.utils.js");
const path = require("path");

const mset_eq_kzg_grandproduct_prover = require("../src/mset_eq_kzg_grandproduct_prover.js");
const mset_eq_kzg_grandproduct_verifier = require("../src/mset_eq_kzg_grandproduct_verifier.js");

describe("Protocols based on grand-products", function () {
    this.timeout(1000000);

    let curve;

    before(async () => {
        curve = await getCurveFromName("bn128");
    });

    after(async () => {
        await curve.terminate();
    });

    it("Should perform the full proving and verifying process of a SNARK for multiset equalities, based on grand-products and KZG", async () => {
        const nBits =  getRandomValue(2, 10);

        const evalsBufferF = getRandomBuffer(2 ** nBits, curve);
        const evalsBufferT = new Uint8Array(evalsBufferF.byteLength);

        evalsBufferT.set(evalsBufferF.slice(0, evalsBufferF.byteLength - curve.Fr.n8), curve.Fr.n8);
        evalsBufferT.set(evalsBufferF.slice(evalsBufferF.byteLength - curve.Fr.n8, evalsBufferF.byteLength), 0);

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await mset_eq_kzg_grandproduct_prover(pTauFilename, evalsBufferF, evalsBufferT);

        const isValid = await mset_eq_kzg_grandproduct_verifier(pTauFilename, proof, nBits);
        assert.ok(isValid);
    });
});