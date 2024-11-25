import { loadPyodide } from "pyodide";

// FIXME: Any type...
function convert(proxy: any) {
  let val = proxy.toJs()
  proxy.destroy()
  return val
}

function global(name: string) {
  let proxy = pyodide.globals.get(name)
  let val = proxy.toJs({dict_converter : Object.fromEntries})
  proxy.destroy()
  return val
}

onmessage = ({ data: { fasta, digestion, csv, missed_cleavages} }) => {
  let motif = "N[^P][TS]"
  let csvFiles = convert(generate_glycopeptides(fasta, digestion, motif, csv, missed_cleavages))
  csvFiles.forEach(([filename, csv]: [any, any]) => {
    const blob = new Blob([csv], { type: "text/csv" });
    const msg  = {
      type: "Result",
      filename,
      blob,
    };
    postMessage(msg);
  });
}

// NOTE: This version needs to match the version of pyodide installed by npm!
// These files can also be hosted locally from `/static` if something ever
// happens to this CDN, but there will be some build-system demons to battle.
const pyodide = await loadPyodide({
  indexURL: "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/",
});

await pyodide.loadPackage(["micropip"]);
const micropip = pyodide.pyimport("micropip");
await micropip.install("theglam==1.0.0");
// If you need to test development version of pgfinder you should build the wheel and copy the resulting .whl to the
// lib/ directory (adajacent to this file), replace the version below and comment out the above (which loads from
// PyPI).
// await micropip.install('./the-glam-0.1.0-py3-none-any.whl');
await pyodide.runPythonAsync("from theglam import *")
const generate_glycopeptides = pyodide.globals.get('generate_glycopeptides');

const msg = {
  type: "Ready",
  // FIXME: Fill this in with a variable that matches the `micropip.install` version!
  version: undefined,
  digestions: global('DIGESTIONS'),
  glycosylation_motifs: global('GLYCOSYLATION_MOTIFS')
};
postMessage(msg);

// console.log(convert(generate_glycopeptides(">A\nPEPTIDE\n>B\nWHOATHERE", "T", "[EI]", "Glycan,Monoisotopic Mass\nABC,300")))
