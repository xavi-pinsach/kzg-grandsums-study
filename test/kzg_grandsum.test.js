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

        let evalsBufferF = new Array(1);
        let evalsBufferT = new Array(1);
        evalsBufferF[0] = getRandomBuffer(2 ** nBits, curve);
        evalsBufferT[0] = new Uint8Array(evalsBufferF[0].byteLength);

        evalsBufferT[0].set(evalsBufferF[0].slice(0, evalsBufferF[0].byteLength - curve.Fr.n8), curve.Fr.n8);
        evalsBufferT[0].set(evalsBufferF[0].slice(evalsBufferF[0].byteLength - curve.Fr.n8, evalsBufferF[0].byteLength), 0);

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT);

        const isValid = await kzg_grandsum_verifier(pTauFilename, proof, nBits);
        assert.ok(isValid);
    });

    it("should perform a Grand Sum Vector ZKG full proving & verifying process", async () => {
        const nPols =  getRandomValue(2, 10);
        const nBits =  getRandomValue(2, 10);

        let evalsBufferF = new Array(nPols);
        let evalsBufferT = new Array(nPols);
        for (let i = 0; i < nPols; i++) {
            evalsBufferF[i] = getRandomBuffer(2 ** nBits, curve);
            evalsBufferT[i] = new Uint8Array(evalsBufferF[i].byteLength);

            evalsBufferT[i].set(evalsBufferF[i].slice(0, evalsBufferF[i].byteLength - curve.Fr.n8), curve.Fr.n8);
            evalsBufferT[i].set(evalsBufferF[i].slice(evalsBufferF[i].byteLength - curve.Fr.n8, evalsBufferF[i].byteLength), 0);
        }

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await kzg_grandsum_prover(pTauFilename, evalsBufferF, evalsBufferT, nPols);

        const isValid = await kzg_grandsum_verifier(pTauFilename, proof, nBits, nPols);
        assert.ok(isValid);
    });
});