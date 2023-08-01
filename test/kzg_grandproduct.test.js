const assert = require("assert");
const { getCurveFromName } = require("ffjavascript");
const {
    getRandomValue,
    getRandomBuffer,
} = require("./test.utils.js");
const path = require("path");

const kzg_grandproduct_prover = require("../src/kzg_grandproduct_prover.js");
const kzg_grandproduct_verifier = require("../src/kzg_grandproduct_verifier.js");

const Logger = require("logplease");
const logger = Logger.create("", { showTimestamp: false });
Logger.setLogLevel("INFO");

describe("grandsums-study", function () {
    this.timeout(1000000);

    let curve;

    before(async () => {
        curve = await getCurveFromName("bn128");
    });

    after(async () => {
        await curve.terminate();
    });

    it("should perform a Grand Product ZKG full proving & verifying process", async () => {
        // const length =  getRandomValue(10);
        const length =  4;

        const evalsBufferA = getRandomBuffer(2 ** length, curve);
        const evalsBufferB = new Uint8Array(evalsBufferA.byteLength);

        // TODO mix the array elements?
        evalsBufferB.set(evalsBufferA.slice(0, evalsBufferA.byteLength - curve.Fr.n8), curve.Fr.n8);
        evalsBufferB.set(evalsBufferA.slice(evalsBufferA.byteLength - curve.Fr.n8, evalsBufferA.byteLength), 0);

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await kzg_grandproduct_prover(evalsBufferA, evalsBufferB, pTauFilename, { logger });

        const isValid = await kzg_grandproduct_verifier(proof, length, pTauFilename, { logger });
        assert.ok(isValid);
    });
});