// Draw phylogenetic tree

let tree;
let treeVStretch = 6;
let treeHStretch = 2;

function addHighlighToTextNode(node) {
  const textElement = $(node);
  let bbox = textElement[0].getBBox();
  let rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  $(rect).attr({
      x: bbox.x - 20,
      y: bbox.y - 10,
      width: bbox.width + 40,
      height: bbox.height + 20,
      fill: "#f35fc8",
  });
  textElement.parent().prepend(rect);
}

function treeIncreaseVStretch(n) {
  treeVStretch += n;
  tree.setVStretch(treeVStretch);
}

function treeIncreaseHStretch(n) {
  treeHStretch += n;
  tree.setHStretch(treeHStretch);
}

document.addEventListener("DOMContentLoaded", function() {
  tree = new TidyTree(newickString, {
    parent: "#phylotree" ,
    layout: 'horizontal',
    mode: 'square',
    vStretch: treeVStretch,
    hStretch: treeHStretch,
  });
  tree.setLeafLabels(true);
  $('#phylotree .tidytree-ruler rect')[0].style.fill = 'none';

  // Highlight query sequence
  tree.eachLeafNode((elem, data) => {
    const label = $(elem).parent().find('text');
    const text = label.text();
    if (text === sampleId) {
      addHighlighToTextNode(label);
    }
    leafName = leafNames[text];
    if (leafName && leafName.species) {
      label.text(leafName.species);
    }
  });
});
