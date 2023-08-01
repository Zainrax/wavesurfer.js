// Spectrogram plugin

import WaveSurfer from 'https://unpkg.com/wavesurfer.js@7/dist/wavesurfer.esm.js'
import Spectrogram from 'https://unpkg.com/wavesurfer.js@7/dist/plugins/spectrogram.esm.js'

// Create an instance of WaveSurfer
const ws = WaveSurfer.create({
  container: '#waveform',
  waveColor: 'rgb(200, 0, 200)',
  progressColor: 'rgb(100, 0, 100)',
  url: '/examples/audio/252697.m4a',
  sampleRate: 48000,
})

// Initialize the Spectrogram plugin
ws.registerPlugin(
  Spectrogram.create({
    labels: true,
    height: 200,
    fftSamples:512,
    labelsColor:    'rgb(100, 0, 100)',
    splitChannels: false,
  }),
)

// Play on click
ws.once('interaction', () => {
  ws.play()
})

/*
<html>
  <div id="waveform"></div>
  <p>
    📖 <a href="https://wavesurfer-js.org/docs/classes/plugins_spectrogram.SpectrogramPlugin">Spectrogram plugin docs</a>
  </p>
</html>
*/
