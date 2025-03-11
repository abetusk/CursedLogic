
import { defineConfig } from 'vite';

// vite.config.js
export default defineConfig({
  optimizeDeps: {
    esbuildOptions: {
      target: 'esnext'
    }
  },
  build: {
    target: 'esnext'
  },
  // more config options ...
})


