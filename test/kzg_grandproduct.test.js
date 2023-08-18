const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const {
    getRandomValue,
    getRandomBuffer,
} = require("./test.utils.js");
const path = require("path");

const kzg_grandproduct_prover = require("../src/kzg_grandproduct_prover.js");
const kzg_grandproduct_verifier = require("../src/kzg_grandproduct_verifier.js");

const logger = require("../logger.js");

describe("grandproduct-study", function () {
    this.timeout(1000000);

    let curve;

    before(async () => {
        curve = await getCurveFromName("bn128");
    });

    after(async () => {
        await curve.terminate();
    });

    it("should perform a Grand Product ZKG full proving & verifying process", async () => {
        // const nBits =  getRandomValue(2, 10);
        const nBits =  getRandomValue(2, 2);

        const evalsBufferA = getRandomBuffer(2 ** nBits, curve);
        const evalsBufferB = new Uint8Array(evalsBufferA.byteLength);

        evalsBufferB.set(evalsBufferA.slice(0, evalsBufferA.byteLength - curve.Fr.n8), curve.Fr.n8);
        evalsBufferB.set(evalsBufferA.slice(evalsBufferA.byteLength - curve.Fr.n8, evalsBufferA.byteLength), 0);

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await kzg_grandproduct_prover(evalsBufferA, evalsBufferB, pTauFilename);

        const isValid = await kzg_grandproduct_verifier(proof, nBits, pTauFilename);
        assert.ok(isValid);
    });
});