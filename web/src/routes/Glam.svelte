<script lang="ts">
	import { FileUpload, type FileUploadApi } from '@skeletonlabs/skeleton-svelte';
	import { onMount } from 'svelte'
	import Glam from '$lib/glam.ts?worker';
	// import { getAcceptAttrString } from "@zag-js/file-utils"
	// Icons
	import IconDropzone from 'lucide-svelte/icons/image-plus';
	import IconFile from 'lucide-svelte/icons/paperclip';
	import IconRemove from 'lucide-svelte/icons/circle-x';
	import fileDownload from 'js-file-download';
	// API Reference
	// FIXME: Pick back up with debugging this... I think there is a place in Zag where
	// the `accept` parameter is getting lost...
	// let apiRef: FileUploadApi | undefined = $state();
	// $inspect(apiRef?.getHiddenInputProps())
	// console.log(getAcceptAttrString({"text/html": [".html", ".htm"], "image/jpeg": [".jpg", ".jpeg"]}))

	// FIXME: Naming?
	let fasta: string | undefined = $state();
	let digestions = $state();
	// FIXME: Is this too hacky?
	let digestion = $state("([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))");
	let csv: string | undefined = $state();

	let digest_settings = $state({
		missed_cleavages:  0,
		min_length:  null,
		max_length:  null,
		semi_enzymatic:  false,
	})

	let glam: Worker;

	onMount(() => {
		glam = new Glam();
		glam.onmessage = ({ data: msg }) => {
			console.log(msg)
			switch (msg.type) {
				case "Ready":
					digestions = msg.digestions;
					break;
				case "Result":
					fileDownload(msg.blob, msg.filename);
					break;
			}
		}
	});

	function runGlam() {
		glam.postMessage({fasta, digestion, csv, missed_cleavages: digest_settings.missed_cleavages});
	}
	
	// FIXME: Obviously replace this `any`
	async function getFile(file: any) {
		return await file.files[0].text();
	}
</script>

<article
	class="card flex w-full max-w-md flex-col border-2 p-4 text-center border-surface-200-800 preset-filled-surface-100-900"
>
	<h1 class="h1">FASTA</h1>
	<FileUpload
		subtext="Upload a protein-containing FASTA file"
		onFileAccept={async (f: any) => (fasta = await getFile(f))}
		classes="w-full"
		accept={{ 'text/x-fasta': ['.fasta', '.fas', '.fa', '.faa', '.mpfa'] }}
	>
		{#snippet iconInterface()}<IconDropzone class="size-8" />{/snippet}
		{#snippet iconFile()}<IconFile class="size-4" />{/snippet}
		{#snippet iconFileRemove()}<IconRemove class="size-4" />{/snippet}
	</FileUpload>

	<form class="flex gap-4 items-center mx-auto w-full max-w-md mt-4">
		<select class="select" bind:value={digestion}>
			<!-- FIXME: Allow selecting multiple enzymes by fusing regex with `|` -->
			{#if digestions === undefined }
				<option>:(</option>
			{:else}
				{#each Object.entries(digestions?) as [name, regex]}
					<option value={regex}>{name}</option>
				{/each}
			{/if}
		</select>
		<input type="number" class="mt-0 input" bind:value={digest_settings.missed_cleavages}/>
	</form>
	
	<!-- FIXME: The `accept` parameter doesn't appear to be working correctly here... -->
	<h1 class="h1">CSV</h1>
	<FileUpload
		subtext="Upload a glycan database file (.csv)"
		onFileAccept={async (f: any) => (csv = await getFile(f))}
		classes="w-full"
		accept="text/csv"
	>
		{#snippet iconInterface()}<IconDropzone class="size-8" />{/snippet}
		{#snippet iconFile()}<IconFile class="size-4" />{/snippet}
		{#snippet iconFileRemove()}<IconRemove class="size-4" />{/snippet}
	</FileUpload>
	<button type="button" class="btn preset-filled mt-2" onclick={runGlam}>Go!</button>
</article>
