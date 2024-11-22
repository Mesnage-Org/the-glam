<script lang="ts">
	import { FileUpload } from '@skeletonlabs/skeleton-svelte';
	// Icons
	import IconDropzone from 'lucide-svelte/icons/image-plus';
	import IconFile from 'lucide-svelte/icons/paperclip';
	import IconRemove from 'lucide-svelte/icons/circle-x';
</script>

<article class="flex flex-col card preset-filled-surface-100-900 border-2 border-surface-200-800 w-full max-w-md p-4 text-center">
	<h1>Card</h1>
	<!-- FIXME: The `accept` parameter doesn't appear to be working correctly here... -->
  <FileUpload
  	accept="image/png"
  	subtext="Upload a glycan database file (.csv)"
  	onFileChange={console.log}
  	onFileReject={console.error}
  	classes="w-full"
  >
  	{#snippet iconInterface()}<IconDropzone class="size-8" />{/snippet}
  	{#snippet iconFile()}<IconFile class="size-4" />{/snippet}
  	{#snippet iconFileRemove()}<IconRemove class="size-4" />{/snippet}
  </FileUpload>
</article>

