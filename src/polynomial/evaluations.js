const { BigBuffer } = require("ffjavascript");

module.exports.Evaluations =  class Evaluations {
    constructor(evaluations, curve, logger) {
        this.eval = evaluations;
        this.curve = curve;
        this.Fr = curve.Fr;
        this.logger = logger;
    }

    static async fromPolynomial(polynomial, extension, curve, logger) {
        const coefficientsN = new BigBuffer(polynomial.length() * extension * curve.Fr.n8);
        coefficientsN.set(polynomial.coef, 0);

        const evaluations = await curve.Fr.fft(coefficientsN);

        return new Evaluations(evaluations, curve, logger);
    }

    getEvaluation(index) {
        const i_n8 = index * this.Fr.n8;

        if (i_n8 + this.Fr.n8 > this.eval.byteLength) {
            throw new Error("Evaluations.getEvaluation() out of bounds");
        }

        return this.eval.slice(i_n8, i_n8 + this.Fr.n8);
    }

    length() {
        let length = this.eval.byteLength / this.Fr.n8;
        if (length !== Math.floor(this.eval.byteLength / this.Fr.n8)) {
            throw new Error("Polynomial evaluations buffer has incorrect size");
        }
        if (0 === length) {
            this.logger.warn("Polynomial has length zero");
        }
        return length;
    }
}