const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const {
    getRandomValue,
    getRandomBuffer,
} = require("./test.utils.js");
const path = require("path");

const kzg_grandsum_prover = require("../src/kzg_grandsum_prover.js");
const kzg_grandsum_verifier = require("../src/kzg_grandsum_verifier.js");

describe("grandsums-study", function () {
    this.timeout(1000000);

    let curve;

    before(async () => {
        curve = await getCurveFromName("bn128");
    });

    after(async () => {
        await curve.terminate();
    });

    it("should perform a Grand Sum ZKG full proving & verifying process", async () => {
        const nBits =  getRandomValue(2, 10);

        const evalsBufferA = getRandomBuffer(2 ** nBits, curve);
        const evalsBufferB = new Uint8Array(evalsBufferA.byteLength);

        evalsBufferB.set(evalsBufferA.slice(0, evalsBufferA.byteLength - curve.Fr.n8), curve.Fr.n8);
        evalsBufferB.set(evalsBufferA.slice(evalsBufferA.byteLength - curve.Fr.n8, evalsBufferA.byteLength), 0);

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await kzg_grandsum_prover(evalsBufferA, evalsBufferB, pTauFilename);

        const isValid = await kzg_grandsum_verifier(proof, nBits, pTauFilename);
        assert.ok(isValid);
    });
});