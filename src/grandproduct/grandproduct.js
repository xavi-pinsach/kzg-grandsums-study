const { BigBuffer } = require("ffjavascript");
const { Polynomial } = require("../polynomial/polynomial");

const logger = require("../../logger.js");

module.exports.ComputeZGrandProductPolynomial = async function ComputeZGrandProductPolynomial(evalsF, evalsT, evalsSelF, evalsSelT, isSelected, challenge, curve) {
    const Fr = curve.Fr;

    logger.info("··· Building the grand-product polynomial Z");

    const n = evalsF.length()

    let numArr = new BigBuffer(evalsF.length() * Fr.n8);
    let denArr = new BigBuffer(evalsF.length() * Fr.n8);

    // Set the first values to 1
    numArr.set(Fr.one, 0);
    denArr.set(Fr.one, 0);

    // Set initial omega
    for (let i = 0; i < n; i++) {
        if ((~i) && (i & 0xFFF === 0)) logger.info(`··· Z evaluation ${i}/${n}`);

        // num := (f + challenge), den := (t + challenge)
        let num = Fr.add(evalsF.getEvaluation(i), challenge);
        let den = Fr.add(evalsT.getEvaluation(i), challenge);

        if (isSelected) {
            // num := self * (num - 1) + 1, den := selt * (den - 1) + 1
            num = Fr.add(Fr.mul(evalsSelF.getEvaluation(i), Fr.sub(num, Fr.one)), Fr.one);
            den = Fr.add(Fr.mul(evalsSelT.getEvaluation(i), Fr.sub(den, Fr.one)), Fr.one);
        }

        numArr.set(num, ((i + 1) % n) * Fr.n8);
        denArr.set(den, ((i + 1) % n) * Fr.n8);
    }

    // Compute the batch inverse of denArr
    denArr = await Fr.batchInverse(denArr);

    // Multiply numArr · denArr where denArr was inverted in the previous command
    let lastVal = Fr.one;
    for (let i = 0; i < n; i++) {
        const i_sFr = ((i + 1) % n) * Fr.n8;

        let s = Fr.mul(numArr.slice(i_sFr, i_sFr + Fr.n8), denArr.slice(i_sFr, i_sFr + Fr.n8));
        lastVal = Fr.mul(s, lastVal);
        numArr.set(lastVal, i_sFr);
    }
    
    // From now on the values saved on numArr will be Z(X) evaluations buffer

    if (!Fr.eq(numArr.slice(0, Fr.n8), Fr.one)) {
        throw new Error("The grand-product polynomial Z is not well calculated");
    }

    // Compute polynomial coefficients z(X) from buffers.Z
    logger.info("··· Computing Z ifft");
    return await Polynomial.fromEvaluations(numArr, curve);
}