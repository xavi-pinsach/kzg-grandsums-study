const { BigBuffer } = require("ffjavascript");
const { Evaluations } = require("./polynomial/evaluations");
const { Polynomial } = require("./polynomial/polynomial");

const logger = require("../logger.js");

module.exports = async function ComputeSGrandSumPolynomial(evaluationsF, evaluationsT, challenge, curve) {
    const evalsF = evaluationsF instanceof Evaluations ? evaluationsF.eval : evaluationsF;
    const evalsT = evaluationsT instanceof Evaluations ? evaluationsT.eval : evaluationsT;

    logger.info("··· Building S Grand Sum polynomial");

    if(evalsF.byteLength !== evalsT.byteLength) {
        throw new Error("Polynomials must have the same size");
    }

    //Check polF and polT buffers length are the same
    const sFr = curve.Fr.n8;
    const n = evalsF.byteLength / sFr;


    let numArr = new BigBuffer(evalsF.byteLength);
    let denArr = new BigBuffer(evalsF.byteLength);

    // Set the first values to 0
    numArr.set(curve.Fr.zero, 0);
    denArr.set(curve.Fr.zero, 0);

    // Set initial omega
    for (let i = 0; i < n; i++) {
        if ((~i) && (i & 0xFFF === 0)) logger.info(`··· S evaluation ${i}/${n}`);
        const i_sFr = i * sFr;

        // f = (f + challenge)
        const f = curve.Fr.add(evalsF.slice(i_sFr, i_sFr + sFr), challenge);

        // t = (t + challenge)
        const t = curve.Fr.add(evalsT.slice(i_sFr, i_sFr + sFr), challenge);

        // 1/f - 1/t = (t - f) / (f * t)
        // num = t - f, den = f * t
        const num = curve.Fr.sub(t, f);
        numArr.set(num, ((i + 1) % n) * sFr);
        const den = curve.Fr.mul(f, t);
        denArr.set(den, ((i + 1) % n) * sFr);
    }

    // Compute the batch inverse of denArr
    denArr = await curve.Fr.batchInverse(denArr);

    // Multiply numArr · denArr where denArr was inverted in the previous command
    let lastVal = curve.Fr.zero;
    for (let i = 0; i < n; i++) {
        const i_sFr = ((i + 1) % n) * sFr;

        let s = curve.Fr.mul(numArr.slice(i_sFr, i_sFr + sFr), denArr.slice(i_sFr, i_sFr + sFr));
        lastVal = curve.Fr.add(s, lastVal);
        numArr.set(lastVal, i_sFr);
    }

    // From now on the values saved on numArr will be S(X) evaluations buffer

    if (!curve.Fr.eq(numArr.slice(0, sFr), curve.Fr.zero)) {
        throw new Error("S(X) grand sum is not well calculated");
    }

    // Compute polynomial coefficients S(X) from buffers.S
    logger.info("··· Computing S ifft");
    const S = await Polynomial.fromEvaluations(numArr, curve, logger);

    delete denArr;

    return S;
}