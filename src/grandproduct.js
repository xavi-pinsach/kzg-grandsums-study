const { BigBuffer } = require("ffjavascript");
const { Evaluations } = require("./polynomial/evaluations");
const { Polynomial } = require("./polynomial/polynomial");

const logger = require("../logger.js");

module.exports = async function ComputeZGrandProductPolynomial(evaluations, challenge, curve) {
    const Fr = curve.Fr;
    const evalsF = evaluations[0][0];
    const evalsT = evaluations[0][1];

    logger.info("··· Building the grand-product polynomial Z");

    if(evalsF.length() !== evalsT.length()) {
        throw new Error("Polynomials must have the same size");
    }

    //Check polF and polT buffers length are the same
    const n = evalsF.length();

    let numArr = new BigBuffer(evalsF.length() * Fr.n8);
    let denArr = new BigBuffer(evalsF.length() * Fr.n8);

    // Set the first values to 1
    numArr.set(Fr.one, 0);
    denArr.set(Fr.one, 0);

    for (let i = 0; i < n; i++) {
        if ((~i) && (i & 0xFFF === 0)) logger.info(`··· Z evaluation ${i}/${n}`);
        const i_sFr = i * Fr.n8;

        // num := (f + challenge), den := (t + challenge)
        let num = Fr.add(evalsF.getEvaluation(i), challenge);
        let den = Fr.add(evalsT.getEvaluation(i), challenge);

        // Multiply the current num or den value with the previous value stored  in numArr or denArr, respectively
        num = Fr.mul(numArr.slice(i_sFr, i_sFr + Fr.n8), num);
        den = Fr.mul(denArr.slice(i_sFr, i_sFr + Fr.n8), den);

        numArr.set(num, ((i + 1) % n) * Fr.n8);
        denArr.set(den, ((i + 1) % n) * Fr.n8);
    }
    // Compute the inverse of denArr to compute in the next step the
    // division numArr/denArr by multiplying num · 1/denArr
    denArr = await Fr.batchInverse(denArr);

    // Perform element-wise multiplication numArr · denArr, where denArr was inverted in the previous command
    for (let i = 0; i < n; i++) {
        const i_sFr = i * Fr.n8;

        const z = Fr.mul(numArr.slice(i_sFr, i_sFr + Fr.n8), denArr.slice(i_sFr, i_sFr + Fr.n8));
        numArr.set(z, i_sFr);
    }
    
    // From now on the values saved on numArr will be Z(X) evaluations buffer

    if (!Fr.eq(numArr.slice(0, Fr.n8), Fr.one)) {
        throw new Error("The grand-product polynomial Z is not well computed");
    }

    // Compute polynomial coefficients z(X) from buffers.Z
    logger.info("··· Computing Z ifft");
    return await Polynomial.fromEvaluations(numArr, curve);
}