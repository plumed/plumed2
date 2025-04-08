var config = {
    startOnLoad:true,
    theme: 'forest',
    securityLevel: 'loose',
    flowchart:{
            useMaxWidth:false,
            htmlLabels:true
        }
};
mermaid.initialize(config);
window.mermaid.init(undefined, document.querySelectorAll('.mermaid'));
