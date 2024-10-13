Actions implemented in PLUMED
-----------------------------
The [actions](actions.md) that can be used within a PLUMED input file are listed below. 

{:#browse-table .display}
| Name | Description |
|:--------:|:--------:|
{% for item in site.data.actionlist %}| [{{ item.name }}]({{ item.path }}) | {{ item.description }} |
{% endfor %}


<script>
$(document).ready(function() {
var table = $('#browse-table').DataTable({
  "dom": '<"search"f><"top"il>rt<"bottom"Bp><"clear">',
  language: { search: '', searchPlaceholder: "Search project..." },
  buttons: [
        'copy', 'excel', 'pdf'
  ],
  "order": [[ 0, "desc" ]]
  });
$('#browse-table-searchbar').keyup(function () {
  table.search( this.value ).draw();
  });
  hu = window.location.search.substring(1);
  searchfor = hu.split("=");
  if( searchfor[0]=="search" ) {
      table.search( searchfor[1] ).draw();
  }
});
</script>
