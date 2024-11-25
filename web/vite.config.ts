import { sveltekit } from '@sveltejs/kit/vite';
import { defineConfig } from 'vite';

export default defineConfig({
	esbuild: {
    supported: {
      'top-level-await': true
    },
  },
  plugins: [sveltekit()],
	worker: {
		format: "es",
	},
});
