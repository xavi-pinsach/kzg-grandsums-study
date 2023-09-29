const { BigBuffer } = require("ffjavascript");

const logger = require("../../logger.js");

class Evaluations {
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

    static fromArray(array, curve) {
        let buffer = new Uint8Array(array.length * curve.Fr.n8);
        for (let i = 0; i < array.length; i++) {
            buffer.set(array[i], i * curve.Fr.n8);
        }
        return new Evaluations(buffer, curve);
    }

    static fromEvals(evals) {
        return new Evaluations(evals.eval.slice(), evals.curve);
    }

    static getOneEvals(length, curve) {
        let buffer = new Uint8Array(length * curve.Fr.n8);
        for (let i = 0; i < length; i++) {
            buffer.set(curve.Fr.one, i * curve.Fr.n8);
        }
        return new Evaluations(buffer, curve);
    }

    static getZeroEvals(length, curve) {
        let buffer = new Uint8Array(length * curve.Fr.n8);
        for (let i = 0; i < length; i++) {
            buffer.set(curve.Fr.zero, i * curve.Fr.n8);
        }
        return new Evaluations(buffer, curve);
    }

    static getRandomEvals(length, curve) {
        let buffer = new Uint8Array(length * curve.Fr.n8);
        for (let i = 0; i < length; i++) {
            buffer.set(curve.Fr.random(), i * curve.Fr.n8);
        }
        return new Evaluations(buffer, curve);
    }

    static getRandomBinEvals(length, curve) {
        let buffer = new Uint8Array(length * curve.Fr.n8);
        for (let i = 0; i < length; i++) {
            const bit = Math.floor(Math.random() * 2);
            buffer.set(bit === 1 ? curve.Fr.one : curve.Fr.zero, i * curve.Fr.n8);
        }
        return new Evaluations(buffer, curve);
    }

    getEvaluation(index) {
        if ((index + 1) * this.Fr.n8 > this.eval.byteLength) {
            throw new Error("Evaluations.getEvaluation() out of bounds");
        }

        return this.eval.slice(index * this.Fr.n8, (index + 1) * this.Fr.n8);
    }

    getEvaluationSequence(start, end) {
        if (start > end) {
            throw new Error("Evaluations.getEvaluationSequence() start index is greater than end index");
        } else if (start === end) {
            throw new Error("Use Evaluations.getEvaluation() instead");
        }

        if (end > this.length() - 1) {
            throw new Error("Evaluations.getEvaluationSequence() end index is out of bounds");
        }

        return this.eval.slice(start * this.Fr.n8, end * this.Fr.n8);
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

    print(name = "f") {
        for (let i = 0; i < this.length(); i++) {
            console.log(`${name}(𝛚^${i}) =`, this.Fr.toString(this.getEvaluation(i)));
        }
    }
}

module.exports = { Evaluations };