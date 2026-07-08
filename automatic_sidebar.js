// for automatic sidebar based on headings
function generateRandomString(length) {
    return Math.random().toString(36).substring(2, 2 + length);
}
/**
 * @param {HTMLElement} el
 */
function make_sidebar(el){
    /** 
     * @type {{title:string , subsections:{ title:string , node:HTMLElement}[] , node:HTMLElement}[]}
     */
    let projects = [];
    for (let child of el.children){
        if (child.hasAttribute("data-no-sidebar")) continue;

        if (child.tagName.toLowerCase() == "h1"){
            projects.push({
                title: child.textContent,
                subsections : [],
                node:child
            });
        }
        if (child.tagName.toLowerCase() == "h2"){
            projects[projects.length-1].subsections.push({
                title: child.textContent,
                node: child
            })
        }
    }
    for (let p of projects){
        let element_id = p.node.getAttribute("id")
        if (element_id == null){
            element_id = generateRandomString(10)
            p.node.setAttribute("id", element_id)
        }

        let projectname = document.createElement("a")
        projectname.innerText = p.title
        projectname.setAttribute("href", "#"+element_id)
        projectname.setAttribute("class", "sidebar-element")
        projectname.classList.add("project-link")
        let sidebar = document.getElementById("sidebar-contents")
        sidebar.appendChild(projectname)

        for (let s of p.subsections){
            let sub_id = s.node.getAttribute("id")
            if (sub_id == null){
                sub_id = generateRandomString(10)
                s.node.setAttribute("id", sub_id)
            }
            let subname = document.createElement("a")
            subname.innerText = s.title
            subname.setAttribute("href", "#"+sub_id)
            subname.setAttribute("class", "sidebar-element")
            subname.classList.add("subsection-link")
            sidebar.appendChild(subname)
        }
    }
}
window.addEventListener("load", ()=>{
make_sidebar(document.getElementById("main-body"))
})