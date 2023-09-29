function getRandomValue(min, max) {
    if(max < 1) return min;

    const random = Math.ceil(Math.random() * max);
    return Math.max(min, random);
}

function getRandomArray(length, curve) {
    let array = [];
    for (let i = 0; i < length; i++) {
        array[i] = curve.Fr.random();
    }
    return array;
}

module.exports = {
    getRandomValue,
    getRandomArray
};