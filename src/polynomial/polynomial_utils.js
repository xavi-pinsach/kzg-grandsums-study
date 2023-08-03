module.exports.computeZHEvaluation = function computeZHEvaluation(curve, x, nBits) {
    const Fr = curve.Fr;

    let xn = x;
    for (let i = 0; i < nBits; i++) {
        xn = Fr.square(xn);
    }

    return Fr.sub(xn, Fr.one);
}

module.exports.computeL1Evaluation = function computeL1Evaluation(curve, x, ZHx, nBits) {
    const Fr = curve.Fr;

    const n = Fr.e(2**nBits);
    const w = Fr.one;

    return Fr.div(Fr.mul(w, ZHx), Fr.mul(n, Fr.sub(x, w)));
}
