const { BigBuffer } = require("ffjavascript");
const { Evaluations } = require("./polynomial/evaluations");
const { Polynomial } = require("./polynomial/polynomial");

module.exports = async function buildZGrandSum(evaluationsF, evaluationsT, challenge, curve, options) {
    const evalsF = evaluationsF instanceof Evaluations ? evaluationsF.eval : evaluationsF;
    const evalsT = evaluationsT instanceof Evaluations ? evaluationsT.eval : evaluationsT;

    const logger = options.logger;

    if (logger) logger.info("··· Building Z Grand Sum polynomial");

    if(evalsF.byteLength !== evalsT.byteLength) {
        throw new Error("Polynomials must have the same size");
    }

    //Check polF and polT buffers length are the same
    const sFr = curve.Fr.n8;
    const n = evalsF.byteLength / sFr;


    let numArr = new BigBuffer(evalsF.byteLength);
    let denArr = new BigBuffer(evalsF.byteLength);

    // Set the first values to 1
    numArr.set(curve.Fr.zero, 0);
    denArr.set(curve.Fr.zero, 1);

    // Set initial omega
    for (let i = 0; i < n; i++) {
        if (logger && (~i) && (i & 0xFFF === 0)) logger.debug(`··· Z evaluation ${i}/${n}`);
        const i_sFr = i * sFr;

        // num := (f + challenge)
        const f = curve.Fr.add(evalsF.slice(i_sFr, i_sFr + sFr), challenge);

        // den := (t + challenge)
        const t = curve.Fr.add(evalsT.slice(i_sFr, i_sFr + sFr), challenge);

        // 1/t - 1/f = (f - t) / (f * t)
        // num = f - t, den = f * t
        const num = curve.Fr.sub(f, t);
        numArr.set(num, ((i + 1) % n) * sFr);
        const den = curve.Fr.mul(f, t);
        denArr.set(den, ((i + 1) % n) * sFr);
    }

    // Compute the batch inverse of denArr
    denArr = await curve.Fr.batchInverse(denArr);

    // Multiply numArr · denArr where denArr was inverted in the previous command
    let lastVal = curve.Fr.zero;
    for (let i = 0; i < n; i++) {
        const i_sFr = i * sFr;

        let z = curve.Fr.mul(numArr.slice(i_sFr, i_sFr + sFr), denArr.slice(i_sFr, i_sFr + sFr));
        lastVal = curve.Fr.add(z, lastVal);
        numArr.set(lastVal, i_sFr);
    }
    
    // From now on the values saved on numArr will be Z(X) evaluations buffer

    if (!curve.Fr.eq(numArr.slice(0, sFr), curve.Fr.zero)) {
        throw new Error("Z(X) grand sum is not well calculated");
    }

    // Compute polynomial coefficients z(X) from buffers.Z
    if (logger) logger.info("··· Computing Z ifft");
    const Z = await Polynomial.fromEvaluations(numArr, curve, logger);

    delete denArr;

    return Z;
}