// const assert = require("assert");
// const { getCurveFromName } = require("ffjavascript");
// const { getRandomPolynomialByLength, getRandomValue } = require("./test.utils.js");
// const { Polynomial } = require("../src/polynomial/polynomial.js");
// const path = require("path");

// const kzg_GS_prover = require("../src/kzg_grandsum_prover.js");
// const kzg_basic_verifier = require("../src/kzg_basic_verifier.js");

// const Logger = require("logplease");
// const logger = Logger.create("", { showTimestamp: false });
// Logger.setLogLevel("INFO");

// describe("grand-sums-study: KZG basic (1 polynomial) test", function () {
//     this.timeout(500000);

//     let curve;

//     before(async () => {
//         curve = await getCurveFromName("bn128");
//     });

//     after(async () => {
//         await curve.terminate();
//     });

//     it.skip("should perform a Grand Product ZKG full proving & verifying process with two polynomials", async () => {
//         // Get a random number of polynomials to be committed between 2 and 5
//         pol = getRandomPolynomialByLength(2, curve.Fr);

//         const pTauFilename = path.join("tmp", "powersOfTau28_hez_final_15.ptau");
//         const proof = await kzg_GS_prover([pol], pTauFilename, { logger });

//         //const isValid = await kzg_basic_verifier(proof, pTauFilename, { logger });
//         assert.ok(isValid);
//     });
// });
