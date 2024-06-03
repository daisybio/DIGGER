// load data from local storage
function getDataFromLocalStorage(key) {
    const jsonData = localStorage.getItem(key);
    return jsonData ? JSON.parse(jsonData) : null;
}

// Function to create an HTML template
function createHtmlTemplate(data) {
    return `
     <div class="row">
      <div class="col-md-12 my-1">
        <div class="card previous-card" onclick="prevAnalysis('${data.value}', '${data.name}')">
          <div class="card-body">
            <h6 class="card-title" style="display: inline; font-weight: bold">${data.name}</h6>
            <p class="card-text mx-1" style="display: inline; color: gray">●</p>
            <p class="card-text" style="display: inline; color: gray">${data.value}</p>
          </div>
        </div>
      </div>
    </div>
    `;
}

// Function to append the template to a div element
function appendTemplateToDiv(template, divId) {
    const container = document.getElementById(divId);
    if (container) {
        container.innerHTML += template;
    } else {
        console.error(`Element with ID ${divId} not found.`);
    }
}

function prevAnalysis(id, name) {
    document.getElementById('previous_analyses_input').value = id;
    document.getElementById('previous_analyses_name').value = name;
    document.getElementById('submit').click();
}