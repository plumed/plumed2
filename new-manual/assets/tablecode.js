$(document).ready(function() {
var table = $('table.display').DataTable({
  "dom": '<"search"f><"top"il>rt<"bottom"Bp><"clear">',
  language: { search: '', searchPlaceholder: "Search table..." },
  buttons: [
        'copy', 'excel', 'pdf'
  ],
  "order": [[ 0, "desc" ]]
  });
$('#browse-table-searchbar').keyup(function () {
  table.search( this.value ).draw();
  });
  var hu = window.location.search.substring(1);
  var searchfor = hu.split("=");
  if( searchfor[0]=="search" ) {
      table.search( searchfor[1] ).draw();
  } else if( searchfor[0]=="tagsearch" ) {
      table.columns(2).search( "\\b" + searchfor[1] + "\\b", true, false, false ).draw();
      displayTagData(searchfor[1]);
  } else if( searchfor[0]=="actions" ) {
      var allact = searchfor[1].split(",");
      var searchstr = "\\b" + allact[0] + "\\b";
      for( let i=1; i<allact.length; i++) {
           searchstr += "|\\b" + allact[i] + "\\b";
      }
      table.columns(0).search( searchstr, true, false, false ).draw();
  }
});

var tag_selected = false;

function displayTagData(name) {
  var valuefield = document.getElementById("tagdisplay");
  if( !tag_selected && name=="no selection" ) {
      valuefield.innerHTML = "";   
      return;
  }
  var tagvalue = document.getElementById(name);
  valuefield.innerHTML = tagvalue.innerHTML;
};

function displayActionWithTags(name) {
  var page = location.href;
  location.replace( page.split("?")[0] + "?tagsearch=" + name);
  tag_selected = true;
};

function clearTagSelection() {
  tag_selected = false;
  var page = location.href;
  location.replace( page.split("?")[0] );
};
