const assert = require("assert");
const path = require("path");
const { getCurveFromName } = require("ffjavascript");
const { getRandomValue } = require("./test.utils.js");
const { Evaluations } = require("../src/polynomial/evaluations.js");

const mset_eq_kzg_grandproduct_prover = require("../src/grandproduct/mset_eq_kzg_prover.js");
const mset_eq_kzg_grandproduct_verifier = require("../src/grandproduct/mset_eq_kzg_verifier.js");

describe("Protocols based on grand-products", function () {
    this.timeout(1000000);

    let curve;
    let Fr;

    before(async () => {
        curve = await getCurveFromName("bn128");
        Fr = curve.Fr;
    });

    after(async () => {
        await curve.terminate();
    });

    it("Should perform the full proving and verifying process of a SNARK for multiset equalities, based on grand-products and KZG", async () => {
        const nBits =  getRandomValue(1, 10);

        const evalsF = Evaluations.getRandomEvals(2 ** nBits, curve);
        const evalsT = Evaluations.fromEvals(evalsF);
        evalsT.setEvaluation(1, evalsF.getEvaluationSequence(0, evalsF.length() - 1));
        evalsT.setEvaluation(0, evalsF.getEvaluation(evalsF.length() - 1));

        const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_11.ptau");
        const proof = await mset_eq_kzg_grandproduct_prover(pTauFilename, evalsF, evalsT);

        const isValid = await mset_eq_kzg_grandproduct_verifier(pTauFilename, proof, nBits);
        assert.ok(isValid);
    });
});