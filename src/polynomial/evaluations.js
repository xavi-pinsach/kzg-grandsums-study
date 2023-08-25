const { BigBuffer } = require("ffjavascript");

const logger = require("../../logger.js");

module.exports.Evaluations =  class Evaluations {
    constructor(evaluations, curve) {
        this.eval = evaluations;
        this.curve = curve;
        this.Fr = curve.Fr;
    }

    static async fromPolynomial(polynomial, extension, curve) {
        const power = Math.ceil(Math.log2(polynomial.length()));
        const length = (1 << power) * extension;
        const coefficientsN = new BigBuffer(length * curve.Fr.n8);
        coefficientsN.set(polynomial.coef, 0);

        const evaluations = await curve.Fr.fft(coefficientsN);

        return new Evaluations(evaluations, curve);
    }

    getEvaluation(index) {
        const i_n8 = index * this.Fr.n8;

        if (i_n8 + this.Fr.n8 > this.eval.byteLength) {
            throw new Error("Evaluations.getEvaluation() out of bounds");
        }

        return this.eval.slice(i_n8, i_n8 + this.Fr.n8);
    }

    setEvaluation(index, value) {
        if (index > this.length() - 1) {
            throw new Error("Evaluation index is out of bounds");
        }

        this.eval.set(value, index * this.Fr.n8);
    }


    length() {
        let length = this.eval.byteLength / this.Fr.n8;
        if (length !== Math.floor(this.eval.byteLength / this.Fr.n8)) {
            throw new Error("Polynomial evaluations buffer has incorrect size");
        }
        if (0 === length) {
            logger.warn("Polynomial has length zero");
        }
        return length;
    }
}