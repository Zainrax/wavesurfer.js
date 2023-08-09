/*!
 * wavesurfer.js spectrogram plugin 5.2.0 (2023-08-09)
 * https://wavesurfer-js.org
 * @license BSD-3-Clause
 */
(function webpackUniversalModuleDefinition(root, factory) {
	if(typeof exports === 'object' && typeof module === 'object')
		module.exports = factory();
	else if(typeof define === 'function' && define.amd)
		define("WaveSurfer", [], factory);
	else if(typeof exports === 'object')
		exports["WaveSurfer"] = factory();
	else
		root["WaveSurfer"] = root["WaveSurfer"] || {}, root["WaveSurfer"]["spectrogram"] = factory();
})(self, () => {
return /******/ (() => { // webpackBootstrap
/******/ 	"use strict";
/******/ 	var __webpack_modules__ = ({

/***/ "./src/plugin/spectrogram/fft.js":
/*!***************************************!*\
  !*** ./src/plugin/spectrogram/fft.js ***!
  \***************************************/
/***/ ((module, exports) => {



Object.defineProperty(exports, "__esModule", ({
  value: true
}));
exports["default"] = FFT;
/* eslint-disable complexity, no-redeclare, no-var, one-var */

/**
 * Calculate FFT - Based on https://github.com/corbanbrook/dsp.js
 *
 * @param {Number} bufferSize Buffer size
 * @param {Number} sampleRate Sample rate
 * @param {Function} windowFunc Window function
 * @param {Number} alpha Alpha channel
 */
function FFT(bufferSize, sampleRate, windowFunc, alpha) {
  this.bufferSize = bufferSize;
  this.sampleRate = sampleRate;
  this.bandwidth = 2 / bufferSize * (sampleRate / 2);
  this.sinTable = new Float32Array(bufferSize);
  this.cosTable = new Float32Array(bufferSize);
  this.windowValues = new Float32Array(bufferSize);
  this.reverseTable = new Uint32Array(bufferSize);
  this.peakBand = 0;
  this.peak = 0;
  var i;
  switch (windowFunc) {
    case 'bartlett':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = 2 / (bufferSize - 1) * ((bufferSize - 1) / 2 - Math.abs(i - (bufferSize - 1) / 2));
      }
      break;
    case 'bartlettHann':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = 0.62 - 0.48 * Math.abs(i / (bufferSize - 1) - 0.5) - 0.38 * Math.cos(Math.PI * 2 * i / (bufferSize - 1));
      }
      break;
    case 'blackman':
      alpha = alpha || 0.16;
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = (1 - alpha) / 2 - 0.5 * Math.cos(Math.PI * 2 * i / (bufferSize - 1)) + alpha / 2 * Math.cos(4 * Math.PI * i / (bufferSize - 1));
      }
      break;
    case 'cosine':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = Math.cos(Math.PI * i / (bufferSize - 1) - Math.PI / 2);
      }
      break;
    case 'gauss':
      alpha = alpha || 0.25;
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = Math.pow(Math.E, -0.5 * Math.pow((i - (bufferSize - 1) / 2) / (alpha * (bufferSize - 1) / 2), 2));
      }
      break;
    case 'hamming':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = 0.54 - 0.46 * Math.cos(Math.PI * 2 * i / (bufferSize - 1));
      }
      break;
    case 'hann':
    case undefined:
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = 0.5 * (1 - Math.cos(Math.PI * 2 * i / (bufferSize - 1)));
      }
      break;
    case 'lanczoz':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = Math.sin(Math.PI * (2 * i / (bufferSize - 1) - 1)) / (Math.PI * (2 * i / (bufferSize - 1) - 1));
      }
      break;
    case 'rectangular':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = 1;
      }
      break;
    case 'triangular':
      for (i = 0; i < bufferSize; i++) {
        this.windowValues[i] = 2 / bufferSize * (bufferSize / 2 - Math.abs(i - (bufferSize - 1) / 2));
      }
      break;
    default:
      throw Error("No such window function '" + windowFunc + "'");
  }
  var limit = 1;
  var bit = bufferSize >> 1;
  var i;
  while (limit < bufferSize) {
    for (i = 0; i < limit; i++) {
      this.reverseTable[i + limit] = this.reverseTable[i] + bit;
    }
    limit = limit << 1;
    bit = bit >> 1;
  }
  for (i = 0; i < bufferSize; i++) {
    this.sinTable[i] = Math.sin(-Math.PI / i);
    this.cosTable[i] = Math.cos(-Math.PI / i);
  }
  this.calculateSpectrum = function (buffer) {
    // Locally scope variables for speed up
    var bufferSize = this.bufferSize,
      cosTable = this.cosTable,
      sinTable = this.sinTable,
      reverseTable = this.reverseTable,
      real = new Float32Array(bufferSize),
      imag = new Float32Array(bufferSize),
      bSi = 2 / this.bufferSize,
      sqrt = Math.sqrt,
      rval,
      ival,
      mag,
      spectrum = new Float32Array(bufferSize / 2);
    var k = Math.floor(Math.log(bufferSize) / Math.LN2);
    if (Math.pow(2, k) !== bufferSize) {
      throw 'Invalid buffer size, must be a power of 2.';
    }
    if (bufferSize !== buffer.length) {
      throw 'Supplied buffer is not the same size as defined FFT. FFT Size: ' + bufferSize + ' Buffer Size: ' + buffer.length;
    }
    var halfSize = 1,
      phaseShiftStepReal,
      phaseShiftStepImag,
      currentPhaseShiftReal,
      currentPhaseShiftImag,
      off,
      tr,
      ti,
      tmpReal;
    for (var i = 0; i < bufferSize; i++) {
      real[i] = buffer[reverseTable[i]] * this.windowValues[reverseTable[i]];
      imag[i] = 0;
    }
    while (halfSize < bufferSize) {
      phaseShiftStepReal = cosTable[halfSize];
      phaseShiftStepImag = sinTable[halfSize];
      currentPhaseShiftReal = 1;
      currentPhaseShiftImag = 0;
      for (var fftStep = 0; fftStep < halfSize; fftStep++) {
        var i = fftStep;
        while (i < bufferSize) {
          off = i + halfSize;
          tr = currentPhaseShiftReal * real[off] - currentPhaseShiftImag * imag[off];
          ti = currentPhaseShiftReal * imag[off] + currentPhaseShiftImag * real[off];
          real[off] = real[i] - tr;
          imag[off] = imag[i] - ti;
          real[i] += tr;
          imag[i] += ti;
          i += halfSize << 1;
        }
        tmpReal = currentPhaseShiftReal;
        currentPhaseShiftReal = tmpReal * phaseShiftStepReal - currentPhaseShiftImag * phaseShiftStepImag;
        currentPhaseShiftImag = tmpReal * phaseShiftStepImag + currentPhaseShiftImag * phaseShiftStepReal;
      }
      halfSize = halfSize << 1;
    }
    for (var i = 0, N = bufferSize / 2; i < N; i++) {
      rval = real[i];
      ival = imag[i];
      mag = bSi * sqrt(rval * rval + ival * ival);
      if (mag > this.peak) {
        this.peakBand = i;
        this.peak = mag;
      }
      spectrum[i] = mag;
    }
    return spectrum;
  };
}
module.exports = exports.default;

/***/ }),

/***/ "./src/plugin/spectrogram/index.js":
/*!*****************************************!*\
  !*** ./src/plugin/spectrogram/index.js ***!
  \*****************************************/
/***/ ((module, exports, __webpack_require__) => {



Object.defineProperty(exports, "__esModule", ({
  value: true
}));
exports["default"] = void 0;
var _fft = _interopRequireDefault(__webpack_require__(/*! ./fft */ "./src/plugin/spectrogram/fft.js"));
function _interopRequireDefault(obj) { return obj && obj.__esModule ? obj : { default: obj }; }
function _typeof(obj) { "@babel/helpers - typeof"; return _typeof = "function" == typeof Symbol && "symbol" == typeof Symbol.iterator ? function (obj) { return typeof obj; } : function (obj) { return obj && "function" == typeof Symbol && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }, _typeof(obj); }
function _regeneratorRuntime() { "use strict"; /*! regenerator-runtime -- Copyright (c) 2014-present, Facebook, Inc. -- license (MIT): https://github.com/facebook/regenerator/blob/main/LICENSE */ _regeneratorRuntime = function _regeneratorRuntime() { return exports; }; var exports = {}, Op = Object.prototype, hasOwn = Op.hasOwnProperty, defineProperty = Object.defineProperty || function (obj, key, desc) { obj[key] = desc.value; }, $Symbol = "function" == typeof Symbol ? Symbol : {}, iteratorSymbol = $Symbol.iterator || "@@iterator", asyncIteratorSymbol = $Symbol.asyncIterator || "@@asyncIterator", toStringTagSymbol = $Symbol.toStringTag || "@@toStringTag"; function define(obj, key, value) { return Object.defineProperty(obj, key, { value: value, enumerable: !0, configurable: !0, writable: !0 }), obj[key]; } try { define({}, ""); } catch (err) { define = function define(obj, key, value) { return obj[key] = value; }; } function wrap(innerFn, outerFn, self, tryLocsList) { var protoGenerator = outerFn && outerFn.prototype instanceof Generator ? outerFn : Generator, generator = Object.create(protoGenerator.prototype), context = new Context(tryLocsList || []); return defineProperty(generator, "_invoke", { value: makeInvokeMethod(innerFn, self, context) }), generator; } function tryCatch(fn, obj, arg) { try { return { type: "normal", arg: fn.call(obj, arg) }; } catch (err) { return { type: "throw", arg: err }; } } exports.wrap = wrap; var ContinueSentinel = {}; function Generator() {} function GeneratorFunction() {} function GeneratorFunctionPrototype() {} var IteratorPrototype = {}; define(IteratorPrototype, iteratorSymbol, function () { return this; }); var getProto = Object.getPrototypeOf, NativeIteratorPrototype = getProto && getProto(getProto(values([]))); NativeIteratorPrototype && NativeIteratorPrototype !== Op && hasOwn.call(NativeIteratorPrototype, iteratorSymbol) && (IteratorPrototype = NativeIteratorPrototype); var Gp = GeneratorFunctionPrototype.prototype = Generator.prototype = Object.create(IteratorPrototype); function defineIteratorMethods(prototype) { ["next", "throw", "return"].forEach(function (method) { define(prototype, method, function (arg) { return this._invoke(method, arg); }); }); } function AsyncIterator(generator, PromiseImpl) { function invoke(method, arg, resolve, reject) { var record = tryCatch(generator[method], generator, arg); if ("throw" !== record.type) { var result = record.arg, value = result.value; return value && "object" == _typeof(value) && hasOwn.call(value, "__await") ? PromiseImpl.resolve(value.__await).then(function (value) { invoke("next", value, resolve, reject); }, function (err) { invoke("throw", err, resolve, reject); }) : PromiseImpl.resolve(value).then(function (unwrapped) { result.value = unwrapped, resolve(result); }, function (error) { return invoke("throw", error, resolve, reject); }); } reject(record.arg); } var previousPromise; defineProperty(this, "_invoke", { value: function value(method, arg) { function callInvokeWithMethodAndArg() { return new PromiseImpl(function (resolve, reject) { invoke(method, arg, resolve, reject); }); } return previousPromise = previousPromise ? previousPromise.then(callInvokeWithMethodAndArg, callInvokeWithMethodAndArg) : callInvokeWithMethodAndArg(); } }); } function makeInvokeMethod(innerFn, self, context) { var state = "suspendedStart"; return function (method, arg) { if ("executing" === state) throw new Error("Generator is already running"); if ("completed" === state) { if ("throw" === method) throw arg; return doneResult(); } for (context.method = method, context.arg = arg;;) { var delegate = context.delegate; if (delegate) { var delegateResult = maybeInvokeDelegate(delegate, context); if (delegateResult) { if (delegateResult === ContinueSentinel) continue; return delegateResult; } } if ("next" === context.method) context.sent = context._sent = context.arg;else if ("throw" === context.method) { if ("suspendedStart" === state) throw state = "completed", context.arg; context.dispatchException(context.arg); } else "return" === context.method && context.abrupt("return", context.arg); state = "executing"; var record = tryCatch(innerFn, self, context); if ("normal" === record.type) { if (state = context.done ? "completed" : "suspendedYield", record.arg === ContinueSentinel) continue; return { value: record.arg, done: context.done }; } "throw" === record.type && (state = "completed", context.method = "throw", context.arg = record.arg); } }; } function maybeInvokeDelegate(delegate, context) { var methodName = context.method, method = delegate.iterator[methodName]; if (undefined === method) return context.delegate = null, "throw" === methodName && delegate.iterator.return && (context.method = "return", context.arg = undefined, maybeInvokeDelegate(delegate, context), "throw" === context.method) || "return" !== methodName && (context.method = "throw", context.arg = new TypeError("The iterator does not provide a '" + methodName + "' method")), ContinueSentinel; var record = tryCatch(method, delegate.iterator, context.arg); if ("throw" === record.type) return context.method = "throw", context.arg = record.arg, context.delegate = null, ContinueSentinel; var info = record.arg; return info ? info.done ? (context[delegate.resultName] = info.value, context.next = delegate.nextLoc, "return" !== context.method && (context.method = "next", context.arg = undefined), context.delegate = null, ContinueSentinel) : info : (context.method = "throw", context.arg = new TypeError("iterator result is not an object"), context.delegate = null, ContinueSentinel); } function pushTryEntry(locs) { var entry = { tryLoc: locs[0] }; 1 in locs && (entry.catchLoc = locs[1]), 2 in locs && (entry.finallyLoc = locs[2], entry.afterLoc = locs[3]), this.tryEntries.push(entry); } function resetTryEntry(entry) { var record = entry.completion || {}; record.type = "normal", delete record.arg, entry.completion = record; } function Context(tryLocsList) { this.tryEntries = [{ tryLoc: "root" }], tryLocsList.forEach(pushTryEntry, this), this.reset(!0); } function values(iterable) { if (iterable) { var iteratorMethod = iterable[iteratorSymbol]; if (iteratorMethod) return iteratorMethod.call(iterable); if ("function" == typeof iterable.next) return iterable; if (!isNaN(iterable.length)) { var i = -1, next = function next() { for (; ++i < iterable.length;) if (hasOwn.call(iterable, i)) return next.value = iterable[i], next.done = !1, next; return next.value = undefined, next.done = !0, next; }; return next.next = next; } } return { next: doneResult }; } function doneResult() { return { value: undefined, done: !0 }; } return GeneratorFunction.prototype = GeneratorFunctionPrototype, defineProperty(Gp, "constructor", { value: GeneratorFunctionPrototype, configurable: !0 }), defineProperty(GeneratorFunctionPrototype, "constructor", { value: GeneratorFunction, configurable: !0 }), GeneratorFunction.displayName = define(GeneratorFunctionPrototype, toStringTagSymbol, "GeneratorFunction"), exports.isGeneratorFunction = function (genFun) { var ctor = "function" == typeof genFun && genFun.constructor; return !!ctor && (ctor === GeneratorFunction || "GeneratorFunction" === (ctor.displayName || ctor.name)); }, exports.mark = function (genFun) { return Object.setPrototypeOf ? Object.setPrototypeOf(genFun, GeneratorFunctionPrototype) : (genFun.__proto__ = GeneratorFunctionPrototype, define(genFun, toStringTagSymbol, "GeneratorFunction")), genFun.prototype = Object.create(Gp), genFun; }, exports.awrap = function (arg) { return { __await: arg }; }, defineIteratorMethods(AsyncIterator.prototype), define(AsyncIterator.prototype, asyncIteratorSymbol, function () { return this; }), exports.AsyncIterator = AsyncIterator, exports.async = function (innerFn, outerFn, self, tryLocsList, PromiseImpl) { void 0 === PromiseImpl && (PromiseImpl = Promise); var iter = new AsyncIterator(wrap(innerFn, outerFn, self, tryLocsList), PromiseImpl); return exports.isGeneratorFunction(outerFn) ? iter : iter.next().then(function (result) { return result.done ? result.value : iter.next(); }); }, defineIteratorMethods(Gp), define(Gp, toStringTagSymbol, "Generator"), define(Gp, iteratorSymbol, function () { return this; }), define(Gp, "toString", function () { return "[object Generator]"; }), exports.keys = function (val) { var object = Object(val), keys = []; for (var key in object) keys.push(key); return keys.reverse(), function next() { for (; keys.length;) { var key = keys.pop(); if (key in object) return next.value = key, next.done = !1, next; } return next.done = !0, next; }; }, exports.values = values, Context.prototype = { constructor: Context, reset: function reset(skipTempReset) { if (this.prev = 0, this.next = 0, this.sent = this._sent = undefined, this.done = !1, this.delegate = null, this.method = "next", this.arg = undefined, this.tryEntries.forEach(resetTryEntry), !skipTempReset) for (var name in this) "t" === name.charAt(0) && hasOwn.call(this, name) && !isNaN(+name.slice(1)) && (this[name] = undefined); }, stop: function stop() { this.done = !0; var rootRecord = this.tryEntries[0].completion; if ("throw" === rootRecord.type) throw rootRecord.arg; return this.rval; }, dispatchException: function dispatchException(exception) { if (this.done) throw exception; var context = this; function handle(loc, caught) { return record.type = "throw", record.arg = exception, context.next = loc, caught && (context.method = "next", context.arg = undefined), !!caught; } for (var i = this.tryEntries.length - 1; i >= 0; --i) { var entry = this.tryEntries[i], record = entry.completion; if ("root" === entry.tryLoc) return handle("end"); if (entry.tryLoc <= this.prev) { var hasCatch = hasOwn.call(entry, "catchLoc"), hasFinally = hasOwn.call(entry, "finallyLoc"); if (hasCatch && hasFinally) { if (this.prev < entry.catchLoc) return handle(entry.catchLoc, !0); if (this.prev < entry.finallyLoc) return handle(entry.finallyLoc); } else if (hasCatch) { if (this.prev < entry.catchLoc) return handle(entry.catchLoc, !0); } else { if (!hasFinally) throw new Error("try statement without catch or finally"); if (this.prev < entry.finallyLoc) return handle(entry.finallyLoc); } } } }, abrupt: function abrupt(type, arg) { for (var i = this.tryEntries.length - 1; i >= 0; --i) { var entry = this.tryEntries[i]; if (entry.tryLoc <= this.prev && hasOwn.call(entry, "finallyLoc") && this.prev < entry.finallyLoc) { var finallyEntry = entry; break; } } finallyEntry && ("break" === type || "continue" === type) && finallyEntry.tryLoc <= arg && arg <= finallyEntry.finallyLoc && (finallyEntry = null); var record = finallyEntry ? finallyEntry.completion : {}; return record.type = type, record.arg = arg, finallyEntry ? (this.method = "next", this.next = finallyEntry.finallyLoc, ContinueSentinel) : this.complete(record); }, complete: function complete(record, afterLoc) { if ("throw" === record.type) throw record.arg; return "break" === record.type || "continue" === record.type ? this.next = record.arg : "return" === record.type ? (this.rval = this.arg = record.arg, this.method = "return", this.next = "end") : "normal" === record.type && afterLoc && (this.next = afterLoc), ContinueSentinel; }, finish: function finish(finallyLoc) { for (var i = this.tryEntries.length - 1; i >= 0; --i) { var entry = this.tryEntries[i]; if (entry.finallyLoc === finallyLoc) return this.complete(entry.completion, entry.afterLoc), resetTryEntry(entry), ContinueSentinel; } }, catch: function _catch(tryLoc) { for (var i = this.tryEntries.length - 1; i >= 0; --i) { var entry = this.tryEntries[i]; if (entry.tryLoc === tryLoc) { var record = entry.completion; if ("throw" === record.type) { var thrown = record.arg; resetTryEntry(entry); } return thrown; } } throw new Error("illegal catch attempt"); }, delegateYield: function delegateYield(iterable, resultName, nextLoc) { return this.delegate = { iterator: values(iterable), resultName: resultName, nextLoc: nextLoc }, "next" === this.method && (this.arg = undefined), ContinueSentinel; } }, exports; }
function asyncGeneratorStep(gen, resolve, reject, _next, _throw, key, arg) { try { var info = gen[key](arg); var value = info.value; } catch (error) { reject(error); return; } if (info.done) { resolve(value); } else { Promise.resolve(value).then(_next, _throw); } }
function _asyncToGenerator(fn) { return function () { var self = this, args = arguments; return new Promise(function (resolve, reject) { var gen = fn.apply(self, args); function _next(value) { asyncGeneratorStep(gen, resolve, reject, _next, _throw, "next", value); } function _throw(err) { asyncGeneratorStep(gen, resolve, reject, _next, _throw, "throw", err); } _next(undefined); }); }; }
function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }
function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, _toPropertyKey(descriptor.key), descriptor); } }
function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); Object.defineProperty(Constructor, "prototype", { writable: false }); return Constructor; }
function _toPropertyKey(arg) { var key = _toPrimitive(arg, "string"); return _typeof(key) === "symbol" ? key : String(key); }
function _toPrimitive(input, hint) { if (_typeof(input) !== "object" || input === null) return input; var prim = input[Symbol.toPrimitive]; if (prim !== undefined) { var res = prim.call(input, hint || "default"); if (_typeof(res) !== "object") return res; throw new TypeError("@@toPrimitive must return a primitive value."); } return (hint === "string" ? String : Number)(input); } /* eslint-enable complexity, no-redeclare, no-var, one-var */
/**
 * @typedef {Object} SpectrogramPluginParams
 * @property {string|HTMLElement} container Selector of element or element in
 * which to render
 * @property {number} fftSamples=512 Number of samples to fetch to FFT. Must be
 * a power of 2.
 * @property {boolean} labels Set to true to display frequency labels.
 * @property {number} noverlap Size of the overlapping window. Must be <
 * fftSamples. Auto deduced from canvas size by default.
 * @property {string} windowFunc='hann' The window function to be used. One of
 * these: `'bartlett'`, `'bartlettHann'`, `'blackman'`, `'cosine'`, `'gauss'`,
 * `'hamming'`, `'hann'`, `'lanczoz'`, `'rectangular'`, `'triangular'`
 * @property {?number} alpha Some window functions have this extra value.
 * (Between 0 and 1)
 * @property {number} pixelRatio=wavesurfer.params.pixelRatio to control the
 * size of the spectrogram in relation with its canvas. 1 = Draw on the whole
 * canvas. 2 = Draw on a quarter (1/2 the length and 1/2 the width)
 * @property {?boolean} deferInit Set to true to manually call
 * `initPlugin('spectrogram')`
 * @property {?number[][]} colorMap A 256 long array of 4-element arrays.
 * Each entry should contain a float between 0 and 1 and specify
 * r, g, b, and alpha.
 */
/**
 * Render a spectrogram visualisation of the audio.
 *
 * @implements {PluginClass}
 * @extends {Observer}
 * @example
 * // es6
 * import SpectrogramPlugin from 'wavesurfer.spectrogram.js';
 *
 * // commonjs
 * var SpectrogramPlugin = require('wavesurfer.spectrogram.js');
 *
 * // if you are using <script> tags
 * var SpectrogramPlugin = window.WaveSurfer.spectrogram;
 *
 * // ... initialising wavesurfer with the plugin
 * var wavesurfer = WaveSurfer.create({
 *   // wavesurfer options ...
 *   plugins: [
 *     SpectrogramPlugin.create({
 *       // plugin options ...
 *     })
 *   ]
 * });
 */
var SpectrogramPlugin = /*#__PURE__*/function () {
  function SpectrogramPlugin(params, ws) {
    var _this = this;
    _classCallCheck(this, SpectrogramPlugin);
    this.drawSpectrogram = /*#__PURE__*/function () {
      var _ref = _asyncToGenerator( /*#__PURE__*/_regeneratorRuntime().mark(function _callee(frequenciesData) {
        var spectrCc, c, pixels, height, width, imageData, i, j, colorMap, redIndex, renderer;
        return _regeneratorRuntime().wrap(function _callee$(_context) {
          while (1) switch (_context.prev = _context.next) {
            case 0:
              if (!isNaN(frequenciesData[0][0])) {
                // data is 1ch [sample, freq] format
                // to [channel, sample, freq] format
                frequenciesData = [frequenciesData];
              }

              // Set the height to fit all channels
              _this.wrapper.style.height = _this.height * frequenciesData.length + 'px';
              _this.canvas.width = _this.width;
              _this.canvas.height = _this.height * frequenciesData.length;
              spectrCc = _this.spectrCc;
              if (spectrCc) {
                _context.next = 7;
                break;
              }
              return _context.abrupt("return");
            case 7:
              c = 0;
            case 8:
              if (!(c < frequenciesData.length)) {
                _context.next = 22;
                break;
              }
              pixels = frequenciesData[c];
              height = pixels[0].length;
              width = pixels.length;
              imageData = new ImageData(width, height);
              for (i = 0; i < width; i++) {
                for (j = 0; j < height; j++) {
                  colorMap = _this.colorMap[pixels[i][j]];
                  redIndex = ((height - 1 - j) * width + i) * 4;
                  imageData.data[redIndex] = colorMap[0] * 255;
                  imageData.data[redIndex + 1] = colorMap[1] * 255;
                  imageData.data[redIndex + 2] = colorMap[2] * 255;
                  imageData.data[redIndex + 3] = colorMap[3] * 255;
                }
              }
              _context.next = 16;
              return createImageBitmap(imageData);
            case 16:
              renderer = _context.sent;
              _context.next = 19;
              return spectrCc.drawImage(renderer, 0, 0, width, height, 0, 0, _this.width, _this.height // destination width, height
              );
            case 19:
              c++;
              _context.next = 8;
              break;
            case 22:
              _this.loadLabels(_this.options.labelsBackground, '12px', '12px', '', _this.options.labelsColor, _this.options.labelsHzColor || _this.options.labelsColor, 'center', '#specLabels', frequenciesData.length);
              _this.emit('ready');
            case 24:
            case "end":
              return _context.stop();
          }
        }, _callee);
      }));
      return function (_x) {
        return _ref.apply(this, arguments);
      };
    }();
    this.params = params;
    this.wavesurfer = ws;
    this.util = ws.util;
    this.frequenciesDataUrl = params.frequenciesDataUrl;
    this._onScroll = function (e) {
      _this.updateScroll(e);
    };
    this._onRender = function () {
      _this.render();
    };
    this._onWrapperClick = function (e) {
      _this._wrapperClickHandler(e);
    };
    this._onReady = function () {
      var drawer = _this.drawer = ws.drawer;
      _this.container = 'string' == typeof params.container ? document.querySelector(params.container) : params.container;
      if (!_this.container) {
        throw Error('No container for WaveSurfer spectrogram');
      }
      if (params.colorMap) {
        if (params.colorMap.length < 256) {
          throw new Error('Colormap must contain 256 elements');
        }
        for (var i = 0; i < params.colorMap.length; i++) {
          var cmEntry = params.colorMap[i];
          if (cmEntry.length !== 4) {
            throw new Error('ColorMap entries must contain 4 values');
          }
        }
        _this.colorMap = params.colorMap;
      } else {
        _this.colorMap = [];
        for (var _i = 0; _i < 256; _i++) {
          var val = (255 - _i) / 256;
          _this.colorMap.push([val, val, val, 1]);
        }
      }
      _this.width = drawer.width;
      _this.pixelRatio = _this.params.pixelRatio || ws.params.pixelRatio;
      _this.fftSamples = _this.params.fftSamples || ws.params.fftSamples || 512;
      _this.height = _this.fftSamples / 2;
      _this.noverlap = params.noverlap;
      _this.windowFunc = params.windowFunc;
      _this.alpha = params.alpha;
      _this.createWrapper();
      _this.createCanvas();
      _this.render();
      drawer.wrapper.addEventListener('scroll', _this._onScroll);
      ws.on('redraw', _this._onRender);
    };
  }
  _createClass(SpectrogramPlugin, [{
    key: "init",
    value: function init() {
      // Check if wavesurfer is ready
      if (this.wavesurfer.isReady) {
        this._onReady();
      } else {
        this.wavesurfer.once('ready', this._onReady);
      }
    }
  }, {
    key: "destroy",
    value: function destroy() {
      this.unAll();
      this.wavesurfer.un('ready', this._onReady);
      this.wavesurfer.un('redraw', this._onRender);
      this.drawer && this.drawer.wrapper.removeEventListener('scroll', this._onScroll);
      this.wavesurfer = null;
      this.util = null;
      this.params = null;
      if (this.wrapper) {
        this.wrapper.removeEventListener('click', this._onWrapperClick);
        this.wrapper.parentNode.removeChild(this.wrapper);
        this.wrapper = null;
      }
    }
  }, {
    key: "createWrapper",
    value: function createWrapper() {
      var prevSpectrogram = this.container.querySelector('spectrogram');
      if (prevSpectrogram) {
        this.container.removeChild(prevSpectrogram);
      }
      var wsParams = this.wavesurfer.params;
      this.wrapper = document.createElement('spectrogram');
      // if labels are active
      if (this.params.labels) {
        var labelsEl = this.labelsEl = document.createElement('canvas');
        labelsEl.classList.add('spec-labels');
        this.drawer.style(labelsEl, {
          left: 0,
          position: 'absolute',
          zIndex: 9,
          height: "".concat(this.height / this.pixelRatio, "px"),
          width: "".concat(55 / this.pixelRatio, "px")
        });
        this.wrapper.appendChild(labelsEl);
        this.loadLabels('rgba(68,68,68,0.5)', '12px', '10px', '', '#fff', '#f7f7f7', 'center', '#specLabels');
      }
      this.drawer.style(this.wrapper, {
        display: 'block',
        position: 'relative',
        userSelect: 'none',
        webkitUserSelect: 'none',
        height: "".concat(this.height / this.pixelRatio, "px")
      });
      if (wsParams.fillParent || wsParams.scrollParent) {
        this.drawer.style(this.wrapper, {
          width: '100%',
          overflowX: 'hidden',
          overflowY: 'hidden'
        });
      }
      this.container.appendChild(this.wrapper);
      this.wrapper.addEventListener('click', this._onWrapperClick);
    }
  }, {
    key: "_wrapperClickHandler",
    value: function _wrapperClickHandler(event) {
      event.preventDefault();
      var relX = 'offsetX' in event ? event.offsetX : event.layerX;
      this.fireEvent('click', relX / this.width || 0);
    }
  }, {
    key: "createCanvas",
    value: function createCanvas() {
      var canvas = this.canvas = this.wrapper.appendChild(document.createElement('canvas'));
      this.spectrCc = canvas.getContext('2d');
      this.util.style(canvas, {
        position: 'absolute',
        zIndex: 4
      });
    }
  }, {
    key: "render",
    value: function render() {
      this.updateCanvasStyle();
      if (this.frequenciesDataUrl) {
        this.loadFrequenciesData(this.frequenciesDataUrl);
      } else {
        this.getFrequencies(this.drawSpectrogram);
      }
    }
  }, {
    key: "updateCanvasStyle",
    value: function updateCanvasStyle() {
      var width = Math.round(this.width / this.pixelRatio) + 'px';
      this.canvas.width = this.width;
      this.canvas.height = this.height;
      this.canvas.style.width = width;
    }
  }, {
    key: "getFrequencies",
    value: function getFrequencies(callback) {
      var fftSamples = this.fftSamples;
      var buffer = this.buffer = this.wavesurfer.backend.buffer;
      var channelOne = buffer.getChannelData(0);
      var bufferLength = buffer.length;
      var sampleRate = buffer.sampleRate;
      var frequencies = [];
      if (!buffer) {
        this.fireEvent('error', 'Web Audio buffer is not available');
        return;
      }
      var noverlap = this.noverlap;
      if (!noverlap) {
        var uniqueSamplesPerPx = buffer.length / this.canvas.width;
        noverlap = Math.max(0, Math.round(fftSamples - uniqueSamplesPerPx));
      }
      var fft = new _fft.default(fftSamples, sampleRate, this.windowFunc, this.alpha);
      var maxSlicesCount = Math.floor(bufferLength / (fftSamples - noverlap));
      var currentOffset = 0;
      while (currentOffset + fftSamples < channelOne.length) {
        var segment = channelOne.slice(currentOffset, currentOffset + fftSamples);
        var spectrum = fft.calculateSpectrum(segment);
        var array = new Uint8Array(fftSamples / 2);
        var j = void 0;
        for (j = 0; j < fftSamples / 2; j++) {
          array[j] = Math.max(-255, Math.log10(spectrum[j]) * 45);
        }
        frequencies.push(array);
        currentOffset += fftSamples - noverlap;
      }
      callback(frequencies, this);
    }
  }, {
    key: "loadFrequenciesData",
    value: function loadFrequenciesData(url) {
      var _this2 = this;
      var request = this.util.fetchFile({
        url: url
      });
      request.on('success', function (data) {
        return _this2.drawSpectrogram(JSON.parse(data), _this2);
      });
      request.on('error', function (e) {
        return _this2.fireEvent('error', e);
      });
      return request;
    }
  }, {
    key: "freqType",
    value: function freqType(freq) {
      return freq >= 1000 ? (freq / 1000).toFixed(1) : Math.round(freq);
    }
  }, {
    key: "unitType",
    value: function unitType(freq) {
      return freq >= 1000 ? 'KHz' : 'Hz';
    }
  }, {
    key: "loadLabels",
    value: function loadLabels(bgFill, fontSizeFreq, fontSizeUnit, fontType, textColorFreq, textColorUnit, textAlign, container) {
      var frequenciesHeight = this.height;
      bgFill = bgFill || 'rgba(68,68,68,0)';
      fontSizeFreq = fontSizeFreq || '12px';
      fontSizeUnit = fontSizeUnit || '10px';
      fontType = fontType || 'Helvetica';
      textColorFreq = textColorFreq || '#fff';
      textColorUnit = textColorUnit || '#fff';
      textAlign = textAlign || 'center';
      container = container || '#specLabels';
      var bgWidth = 55;
      var getMaxY = frequenciesHeight || 512;
      var labelIndex = 5 * (getMaxY / 256);
      var freqStart = 0;
      var step = (this.wavesurfer.backend.ac.sampleRate / 2 - freqStart) / labelIndex;

      // prepare canvas element for labels
      var ctx = this.labelsEl.getContext('2d');
      this.labelsEl.height = this.height;
      this.labelsEl.width = bgWidth;
      var scale = this.height / (this.wavesurfer.backend.ac.sampleRate / 2);
      if (ctx) {
        // fill background
        ctx.fillStyle = bgFill;
        ctx.fillRect(0, 0, bgWidth, getMaxY);
        ctx.fill();
        var i;

        // render labels
        for (i = 0; i <= labelIndex; i++) {
          ctx.textAlign = textAlign;
          ctx.textBaseline = 'middle';
          var freq = freqStart + step * i;
          var label = this.freqType(freq);
          var units = this.unitType(freq);
          var yLabelOffset = 2;
          var x = 16;
          var y = freq * scale;
          y = this.height - y - 10;
          ctx.fillStyle = textColorUnit;
          ctx.font = fontSizeUnit + ' ' + fontType;
          ctx.fillText(units, x + 24, y);
          // freq label
          ctx.fillStyle = textColorFreq;
          ctx.font = fontSizeFreq + ' ' + fontType;
          ctx.fillText(label, x, y);
        }
      }
    }
  }, {
    key: "updateScroll",
    value: function updateScroll(e) {
      if (this.wrapper) {
        this.wrapper.scrollLeft = e.target.scrollLeft;
      }
    }
  }, {
    key: "resample",
    value: function resample(oldMatrix) {
      var columnsNumber = this.width;
      var newMatrix = [];
      var oldPiece = 1 / oldMatrix.length;
      var newPiece = 1 / columnsNumber;
      var i;
      for (i = 0; i < columnsNumber; i++) {
        var column = new Array(oldMatrix[0].length);
        var j = void 0;
        for (j = 0; j < oldMatrix.length; j++) {
          var oldStart = j * oldPiece;
          var oldEnd = oldStart + oldPiece;
          var newStart = i * newPiece;
          var newEnd = newStart + newPiece;
          var overlap = oldEnd <= newStart || newEnd <= oldStart ? 0 : Math.min(Math.max(oldEnd, newStart), Math.max(newEnd, oldStart)) - Math.max(Math.min(oldEnd, newStart), Math.min(newEnd, oldStart));
          var k = void 0;
          /* eslint-disable max-depth */
          if (overlap > 0) {
            for (k = 0; k < oldMatrix[0].length; k++) {
              if (column[k] == null) {
                column[k] = 0;
              }
              column[k] += overlap / newPiece * oldMatrix[j][k];
            }
          }
          /* eslint-enable max-depth */
        }

        var intColumn = new Uint8Array(oldMatrix[0].length);
        var m = void 0;
        for (m = 0; m < oldMatrix[0].length; m++) {
          intColumn[m] = column[m];
        }
        newMatrix.push(intColumn);
      }
      return newMatrix;
    }
  }], [{
    key: "create",
    value:
    /**
     * Spectrogram plugin definition factory
     *
     * This function must be used to create a plugin definition which can be
     * used by wavesurfer to correctly instantiate the plugin.
     *
     * @param  {SpectrogramPluginParams} params Parameters used to initialise the plugin
     * @return {PluginDefinition} An object representing the plugin.
     */
    function create(params) {
      return {
        name: 'spectrogram',
        deferInit: params && params.deferInit ? params.deferInit : false,
        params: params,
        staticProps: {
          FFT: _fft.default
        },
        instance: SpectrogramPlugin
      };
    }
  }]);
  return SpectrogramPlugin;
}();
exports["default"] = SpectrogramPlugin;
module.exports = exports.default;

/***/ })

/******/ 	});
/************************************************************************/
/******/ 	// The module cache
/******/ 	var __webpack_module_cache__ = {};
/******/ 	
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/ 		// Check if module is in cache
/******/ 		var cachedModule = __webpack_module_cache__[moduleId];
/******/ 		if (cachedModule !== undefined) {
/******/ 			return cachedModule.exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = __webpack_module_cache__[moduleId] = {
/******/ 			// no module.id needed
/******/ 			// no module.loaded needed
/******/ 			exports: {}
/******/ 		};
/******/ 	
/******/ 		// Execute the module function
/******/ 		__webpack_modules__[moduleId](module, module.exports, __webpack_require__);
/******/ 	
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/ 	
/************************************************************************/
/******/ 	
/******/ 	// startup
/******/ 	// Load entry module and return exports
/******/ 	// This entry module is referenced by other modules so it can't be inlined
/******/ 	var __webpack_exports__ = __webpack_require__("./src/plugin/spectrogram/index.js");
/******/ 	
/******/ 	return __webpack_exports__;
/******/ })()
;
});
//# sourceMappingURL=wavesurfer.spectrogram.js.map