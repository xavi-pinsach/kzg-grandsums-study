const { BigBuffer } = require("ffjavascript");
const { Evaluations } = require("./polynomial/evaluations");
const { Polynomial } = require("./polynomial/polynomial");

const logger = require("../logger.js");

module.exports = async function ComputeZGrandProductPolynomial(evaluationsF, evaluationsT, challenge, curve) {
    const evalsF = evaluationsF instanceof Evaluations ? evaluationsF.eval : evaluationsF;
    const evalsT = evaluationsT instanceof Evaluations ? evaluationsT.eval : evaluationsT;

    logger.info("··· Building the grand-roduct polynomial Z");

    if(evalsF.byteLength !== evalsT.byteLength) {
        throw new Error("Polynomials must have the same size");
    }

    //Check polF and polT buffers length are the same
    const sFr = curve.Fr.n8;
    const n = evalsF.byteLength / sFr;

    let numArr = new BigBuffer(evalsF.byteLength);
    let denArr = new BigBuffer(evalsF.byteLength);

    // Set the first values to 1
    numArr.set(curve.Fr.one, 0);
    denArr.set(curve.Fr.one, 0);

    for (let i = 0; i < n; i++) {
        if ((~i) && (i & 0xFFF === 0)) logger.info(`··· Z evaluation ${i}/${n}`);
        const i_sFr = i * sFr;

        // num := (f + challenge)
        let num = curve.Fr.add(evalsF.slice(i_sFr, i_sFr + sFr), challenge);

        // den := (t + challenge)
        let den = curve.Fr.add(evalsT.slice(i_sFr, i_sFr + sFr), challenge);

        // Multiply current num value with the previous one saved in numArr
        num = curve.Fr.mul(numArr.slice(i_sFr, i_sFr + sFr), num);
        numArr.set(num, ((i + 1) % n) * sFr);

        // Multiply current den value with the previous one saved in denArr
        den = curve.Fr.mul(denArr.slice(i_sFr, i_sFr + sFr), den);
        denArr.set(den, ((i + 1) % n) * sFr);
    }
    // Compute the inverse of denArr to compute in the next step the
    // division numArr/denArr by multiplying num · 1/denArr
    denArr = await curve.Fr.batchInverse(denArr);

    // Multiply numArr · denArr where denArr was inverted in the previous command
    for (let i = 0; i < n; i++) {
        const i_sFr = i * sFr;

        const z = curve.Fr.mul(numArr.slice(i_sFr, i_sFr + sFr), denArr.slice(i_sFr, i_sFr + sFr));
        numArr.set(z, i_sFr);
    }
    
    // From now on the values saved on numArr will be Z(X) evaluations buffer

    if (!curve.Fr.eq(numArr.slice(0, sFr), curve.Fr.one)) {
        throw new Error("The grand-product Z is not well computed");
    }

    // Compute polynomial coefficients z(X) from buffers.Z
    logger.info("··· Computing Z ifft");
    const Z = await Polynomial.fromEvaluations(numArr, curve, logger);

    delete denArr;

    return Z;
}