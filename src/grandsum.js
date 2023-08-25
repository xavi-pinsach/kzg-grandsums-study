const { BigBuffer } = require("ffjavascript");
const { Polynomial } = require("./polynomial/polynomial");

const logger = require("../logger.js");

module.exports = async function ComputeSGrandSumPolynomial(evalsF, evalsT, challenge, curve) {
    const Fr = curve.Fr;

    logger.info("··· Building the grand-sum polynomial S");

    if(evalsF.length() !== evalsT.length()) {
        throw new Error("Evaluations must have the same size");
    }

    //Check polF and polT buffers length are the same
    const n = evalsF.length()

    let numArr = new BigBuffer(evalsF.length() * Fr.n8);
    let denArr = new BigBuffer(evalsF.length() * Fr.n8);

    // Set the first values to 0
    numArr.set(Fr.zero, 0);
    denArr.set(Fr.zero, 0);

    // Set initial omega
    for (let i = 0; i < n; i++) {
        if ((~i) && (i & 0xFFF === 0)) logger.info(`··· S evaluation ${i}/${n}`);
        const i_sFr = i * Fr.n8;

        // f = (f + challenge), t = (t + challenge)
        const f = Fr.add(evalsF.getEvaluation(i), challenge);
        const t = Fr.add(evalsT.getEvaluation(i), challenge);

        // 1/f - 1/t = (t - f) / (f * t)
        // num = t - f, den = f * t
        const num = Fr.sub(t, f);
        const den = Fr.mul(f, t);

        numArr.set(num, ((i + 1) % n) * Fr.n8);
        denArr.set(den, ((i + 1) % n) * Fr.n8);
    }

    // Compute the batch inverse of denArr
    denArr = await Fr.batchInverse(denArr);

    // Multiply numArr · denArr where denArr was inverted in the previous command
    let lastVal = Fr.zero;
    for (let i = 0; i < n; i++) {
        const i_sFr = ((i + 1) % n) * Fr.n8;

        let s = Fr.mul(numArr.slice(i_sFr, i_sFr + Fr.n8), denArr.slice(i_sFr, i_sFr + Fr.n8));
        lastVal = Fr.add(s, lastVal);
        numArr.set(lastVal, i_sFr);
    }

    // From now on the values saved on numArr will be S(X) evaluations buffer

    if (!Fr.eq(numArr.slice(0, Fr.n8), Fr.zero)) {
        throw new Error("The grand-sum polynomial S is not well calculated");
    }

    // Compute polynomial coefficients S(X) from buffers.S
    logger.info("··· Computing S ifft");
    return await Polynomial.fromEvaluations(numArr, curve, logger);
}