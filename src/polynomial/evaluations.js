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

    static allOnes(length, curve) {
        let oneBuffer = new Uint8Array(length * curve.Fr.n8);
        for (let i = 0; i < length; i++) {
            oneBuffer.set(curve.Fr.one, i * curve.Fr.n8);
        }
        return oneBuffer;
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

    isEqual(other) {
        if (this.length() !== other.length()) {
            return false;
        }
        const result = Buffer.compare(this.eval, other.eval);
        return result === 0 ? true : false;
    }

    isAllZeros() {
        const zeroBuffer = new Uint8Array(this.length() * this.Fr.n8);
        return this.isEqual(new Evaluations(zeroBuffer, this.curve));
    }

    isAllOnes() {
        let oneBuffer = new Uint8Array(this.length() * this.Fr.n8);
        for (let i = 0; i < this.length(); i++) {
            oneBuffer.set(this.Fr.one, i * this.Fr.n8);
        }
        return this.isEqual(new Evaluations(oneBuffer, this.curve));
    }
}