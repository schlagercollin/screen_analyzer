var global_dataframes = "Placeholder";
var display_type = "Placeholder";
var progress_counter = 0;

var analysis = "Placeholder";

function loadJSON (jsonFile) {
    var result = {};
    $.getJSON(jsonFile, function (data) {
        $.each(data, function (index, value) {
           result[index] = value;
        });
    });
    return result;
};

// var graph_config = loadJSON("static/js/graph_config.json");

var graph_config = {
    "Top Sorted Gene Counts": {
        "data": {
            "data_name": "Counts",
            "data_level": "Gene"
        },
        "x": {
            "title": "Target Gene Symbol",
            "column_name": "Target Gene Symbol"
        },
        "y": {
            "title": "Rep1: Top Sorted Gene Population",
            "column_name": "Rep1: Top Sorted Population"
        },
        "text_data": "Target Gene Symbol",
        "sort_by": {
            "value": "Rep1: Top Sorted Population",
            "ascending": false
        }
    },
    "Bottom Sorted Gene Counts": {
        "data": {
            "data_name": "Counts",
            "data_level": "Gene"
        },
        "x": {
            "title": "Target Gene Symbol",
            "column_name": "Target Gene Symbol"
        },
        "y": {
            "title": "Rep1: Bottom Sorted Gene Population",
            "column_name": "Rep1: Bottom Sorted Population"
        },
        "text_data": "Target Gene Symbol",
        "sort_by": {
            "value": "Rep1: Bottom Sorted Population",
            "ascending": false
        }
    },
    "Unsorted Gene Counts": {
        "data": {
            "data_name": "Counts",
            "data_level": "Gene"
        },
        "x": {
            "title": "Target Gene Symbol",
            "column_name": "Target Gene Symbol"
        },
        "y": {
            "title": "Rep1: Unsorted Gene Population",
            "column_name": "Rep1: Unsorted Population"
        },
        "text_data": "Target Gene Symbol",
        "sort_by": {
            "value": "Rep1: Unsorted Population",
            "ascending": false
        }
    },
    "Mageck Top Results": {
        "data": {
            "data_name": "Mageck Top",
            "data_level": "Gene"
        },
        "x": {
            "title": "pos|lfc",
            "column_name": "pos|lfc"
        },
        "y": {
            "title": "pos|-log(p-value)",
            "column_name": "pos|-log(p-value)"
        },
        "text_data": "id",
        "sort_by": {
            "value": "pos|-log(p-value)",
            "ascending": false
        }
    },
    "Mageck Bottom Results": {
        "data": {
            "data_name": "Mageck Bottom",
            "data_level": "Gene"
        },
        "x": {
            "title": "pos|lfc",
            "column_name": "pos|lfc"
        },
        "y": {
            "title": "pos|-log(p-value)",
            "column_name": "pos|-log(p-value)"
        },
        "text_data": "id",
        "sort_by": {
            "value": "pos|-log(p-value)",
            "ascending": false
        }
    },
    "Ratio Genes": {
        "data": {
            "data_name": "Ratio",
            "data_level": "Gene"
        },
        "x": {
            "title": "mean log2(ratio)",
            "column_name": "Mean log2(ratio)"
        },
        "y": {
            "title": "mean zscore",
            "column_name": "Mean ZScore"
        },
        "text_data": "Target Gene Symbol",
        "sort_by": {
            "value": "Mean ZScore",
            "ascending": false
        }
    },
    "Ratio Guides": {
        "data": {
            "data_name": "Ratio",
            "data_level": "Guides"
        },
        "x": {
            "title": "log2( mean abundance )",
            "column_name": "log_MA"
        },
        "y": {
            "title": "ZScore",
            "column_name": "ZScore"
        },
        "text_data": "Target Gene Symbol",
        "sort_by": {
            "value": "ZScore",
            "ascending": false
        }
    }
}

String.prototype.capitalize = function(){
       return this.replace( /(^|\s)([a-z])/g , function(m,p1,p2){ return p1+p2.toUpperCase(); } );
};

$.fn.dataTable.render.ellipsis = function ( cutoff, wordbreak, escapeHtml ) {
    var esc = function ( t ) {
        return t
            .replace( /&/g, '&amp;' )
            .replace( /</g, '&lt;' )
            .replace( />/g, '&gt;' )
            .replace( /"/g, '&quot;' );
    };

    return function ( d, type, row ) {
        // Order, search and type get the original data
        if ( type !== 'display' ) {
            return d;
        }

        if ( typeof d !== 'number' && typeof d !== 'string' ) {
            return d;
        }

        d = d.toString(); // cast numbers

        if ( d.length < cutoff ) {
            return d;
        }

        var shortened = d.substr(0, cutoff-1);

        // Find the last white space character in the string
        if ( wordbreak ) {
            shortened = shortened.replace(/\s([^\s]*)$/, '');
        }

        // Protect against uncontrolled HTML input
        if ( escapeHtml ) {
            shortened = esc( shortened );
        }

        return  '<span class="full_length_str">'+esc(d)+'</span>' +
                '<span class="shortened ellipsis" title="'+esc(d)+'">'+shortened+'&#8230;</span>';
    };
};


function load_analysis(formData, renderPlots){
    // formData needs "analysis_name" key
    console.log("Render Plots? ", renderPlots);
    $.ajax({
        url: "/analysis/load",
        type: 'POST',
        data: formData,
        success: function(result) {
            $("button#load_result").html("<i class='fas fa-chart-line'></i> Load Analysis");
            analyses = result.analysis_info["Analyses Queued"];
            if (renderPlots) {
                if (analyses["Ratio"]){
                    generatePlot(result.name, "Ratio Guides");
                };
                if (analyses["Mageck"]){
                    generatePlot(result.name, "Mageck Top Results");
                    generatePlot(result.name, "Mageck Bottom Results");
                };
            };
            createDataTable(result.data, result.columns);
            $("#collapseOne").collapse('toggle');
        },
        error: function(error) {
            // From responseText (HTML string), extract error message from title element
            var error_HTML_response = $(error.responseText);
            console.log(error_HTML_response);
            var error = $(error_HTML_response).filter("title").html();
            alert("Error loading files from analysis directory. " +
            "Server says: "+ error);
            $("button#load_result").html("<i class='fas fa-chart-line'></i> Load Analysis");
        },
        processData: false,
        contentType: false
    });
};

function load_specific_data(analysis_name, analysis_type, level, columns, sort_by, ascending, datapoints=1000){
    console.log(columns);
    var formData = new FormData();
    formData.set("analysis_name", analysis_name);
    formData.set("analysis_type", analysis_type);
    formData.set("analysis_level", level);
    formData.set("analysis_columns", columns);
    formData.set("analysis_sort_by", sort_by);
    formData.set("analysis_ascending", ascending);
    formData.set("analysis_datapoints", datapoints);
    var promise = Promise.resolve(
                $.ajax({
                    url: "/analysis/load_specific",
                    type: 'POST',
                    data: formData,
                    success: function(result) {
                        result = result.data;
                    },
                    error: function(error) {
                        // From responseText (HTML string), extract error message from title element
                        var error_HTML_response = $(error.responseText);
                        console.log(error_HTML_response);
                        var error = $(error_HTML_response).filter("title").html();
                        alert("Error loading files from analysis directory. " +
                        "Server says: "+ error);
                    },
                    processData: false,
                    contentType: false
                }));
    return promise
};



function generatePlot(analysis_name, plot_type){
    console.log("Generating plot for", plot_type);
    var div_name = plot_type.replace(/ /g,"_");
    var plot_config = graph_config[plot_type];
    var title = plot_type;
    var data_name = plot_config["data"]["data_name"];
    var data_level = plot_config["data"]["data_level"];
    var sort_by = plot_config["sort_by"]["value"];
    var ascending = plot_config["sort_by"]["ascending"];
    var x_data_promise = load_specific_data(analysis_name, data_name, data_level, plot_config["x"]["column_name"], sort_by, ascending);
    var y_data_promise = load_specific_data(analysis_name, data_name, data_level, plot_config["y"]["column_name"], sort_by, ascending);
    var text_data_promise = load_specific_data(analysis_name, data_name, data_level, plot_config["text_data"], sort_by, ascending);
    var xaxis_title = plot_config["x"]["title"];
    var yaxis_title = plot_config["y"]["title"];
    Promise.all([x_data_promise, y_data_promise, text_data_promise]).then(function(values) {
        var x_data = values[0]["data"];
        var y_data = values[1]["data"];
        var text_data = values[2]["data"];
        // console.log(y_data);
        createPlot(div_name, x_data, y_data, text_data, title, xaxis_title, yaxis_title);
    });

}

function createPlot(div_name, x_data, y_data, text_data, title, xaxis_title, yaxis_title){
    div_name += "_plot";
    if ($("#"+div_name).length){
        Plotly.purge(div_name);
    } else {
        $("#plotsDiv").append("<div id='"+div_name+"'></div>");
        $("#plotsDiv").addClass("myPlotContainer")
    }
    var trace = {
        x: x_data,
        y: y_data,
        mode: 'markers',
        text: text_data,
        marker: {
            line: {width: 1}
        }
    };
    var data = [ trace ];
    var layout = {
        autosize: true,
        //width: 1300,
        //height: 600,
        title: title,
        xaxis: {
            autorange: true,
            title: xaxis_title
        },
        yaxis: {
            //type: 'log',
            autorange: true,
            title: yaxis_title
        },
        hovermode:'closest'
    }
    Plotly.newPlot(div_name, data, layout); //create plot
}

function format_columns(columns) {
    formatted = [];
    ordered = [
        'Target Gene Symbol',
        'Description'
    ]
    column_types = {

    }
    for (val of ordered) {
        formatted.push({"mData": val, "sTitle": val, "sClass": "main_data", "bVisible": true});
    }
    var remaining = $(columns).not(ordered).get();
    for (val of remaining) {
        if (val.startsWith("Rep")){
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "count_data"});
        } else if (val.startsWith("pos|")){
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "mageck_data pos_mageck", "sType": "numeric"});
        } else if (val.startsWith("neg|")){
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "mageck_data neg_mageck", "sType": "numeric"});
        } else if (val.startsWith("Mean")){
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "ratio_data", "sType": "numeric"});
        } else {
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "misc_data"});
        }
    }
    return formatted;
}

function createDataTable(data, columns){
    // Destroy table if initialized
    if ( $.fn.DataTable.isDataTable('#resultTable') ) {
      $('#resultTable').DataTable().destroy();
      $('#resultTable tbody').empty();
    }

    // Format columns
    var formatted_columns = format_columns(columns);
    $("#resultTable").removeClass("invisible");
    // Initialize table
    var table = $('#resultTable').DataTable( {
        aaData: data,
        aoColumns: formatted_columns,
        lengthChange: false,
        colReorder: true,
        lengthMenu: [ 5, 10, 25, 50, 75, 100, 200],
        pageLength: 10,
        buttons: [
            {
                extend: 'collection',
                text: 'Display',
                dropup: true,
                buttons: [
                            {
                                extend: 'colvisGroup',
                                text: 'Default',
                                show: ['.main_data'],
                                hide: ['.misc_data', '.count_data', '.pos_mageck', '.neg_mageck', '.ratio_data']
                            },
                            {
                                extend: 'colvisGroup',
                                text: 'Count Data',
                                show: ['.main_data', '.count_data'],
                                hide: ['.misc_data', '.neg_mageck', '.pos_mageck']
                            },
                            {
                                extend: 'colvisGroup',
                                text: 'Mageck | Pos',
                                show: ['.main_data', '.pos_mageck'],
                                hide: ['.misc_data', '.count_data', '.neg_mageck']
                            },
                            {
                                extend: 'colvisGroup',
                                text: 'Mageck | Neg',
                                show: ['.main_data', '.neg_mageck'],
                                hide: ['.misc_data', '.count_data', '.pos_mageck']
                            },
                            {
                                extend: 'colvisGroup',
                                text: 'Ratio',
                                show: ['.main_data','.ratio_data'],
                                hide: ['.misc_data', '.count_data', '.pos_mageck', '.neg_mageck']
                            },
                            {
                                extend: 'colvis',
                                text: 'Custom'
                            }
                        ]
            },
            'pageLength'
        ],
        columnDefs: [ {
            targets: "_all",
            render: $.fn.dataTable.render.ellipsis( 50, true )
        }]
    } );

    table.buttons().container()
        .appendTo( $('div.eight.column:eq(0)', table.table().container()) );
}

function createExpandedRowMarkup ( d ) {
    // `d` is the original data object for the row
    left_markup = "";
    center_markup = "";
    right_markup = "";
    var i = 0;
    for ( const [ key, value ] of Object.entries(d)) {
        markup = `
            <tr>
                <td style="font-style: bold"> ${ key } </td>
                <td> ${ value } </td>
            </tr>
        `
        switch (i%3) {
            case 0:
                left_markup+=markup;
                break;
            case 1:
                center_markup+=markup;
                break;
            case 2:
                right_markup+=markup;
                break;
        }
        i += 1;
    };
    var gene_name = d["Target Gene Symbol"]

    full_markup = `
    <div class="slider">
        <table cellpadding="5" cellspacing="0" border="0" class="detailTable">
            ${ left_markup }
        </table>
        <table cellpadding="5" cellspacing="0" border="0" class="detailTable">
            ${ center_markup }
        </table>
        <table cellpadding="5" cellspacing="0" border="0" class="detailTable">
            ${ right_markup }
        </table>
    </div>
    `
    return full_markup;
}

$( document ).ready( function() {


    class Analysis {
        constructor(analysis_config) {
            this.analysis_config = analysis_config;
            this.info = {
                "Name": analysis_config["Metadata"]["Analysis Name"],
                "Timestamp": analysis_config["Metadata"]["Timestamp"],
                "Replicates": analysis_config["Metadata"]["Replicates"],
                "Analyses Performed": analysis_config["Analyses Queued"]
            };
            var time_obj = new Date(parseInt(this.info["Timestamp"]));
            this.info["Timestamp"] = time_obj.toLocaleString();
        }
    }

    $("select.file_select").change(function(){
        var file_upload_field = $(this).siblings("input.file_upload");
        if ($(this).val() == "Upload your own") {
            $(this).animate({
                'max-width': 17
            }, 200, function() {
                console.log("Test");
                file_upload_field.fadeIn(200);
            });
        } else {
            file_upload_field.fadeOut(0);
            $(this).animate({
                'max-width': 600
            }, 200, function () {

            })
        }
    });

    $("form#load_submit").submit(function(e) {
        //Submit load data request (controller.analysis_load())
        e.preventDefault();
        // Destroy table if initialized
        $("#resultTable").addClass("invisible");
        if ( $.fn.DataTable.isDataTable('#resultTable') ) {
          $('#resultTable').DataTable().destroy();
          $('#resultTable tbody').empty();
        }
        $("button#load_result").text("Loading...");
        var formData = new FormData(this);
        var renderPlots = $("#renderPlotsBool").is(":checked");
        load_analysis(formData, renderPlots);
    });

    //
    $(document.body).on('dblclick', '#resultTable tbody tr[role="row"]', function () {
        var table = $("#resultTable").DataTable();
        var tr = $(this);
        var row = table.row( tr );
        var data = row.data();
        var toggleCell = $(this).find(".toggle-details-row");
        // alert( 'You clicked on '+data['Target Gene Symbol']+'\'s row' );
        //
        // console.log("Modal should pop up.");
        // $('.modal .modal-body').html(createExpandedRowMarkup(row.data()));
        // $('.modal').modal('toggle');

        if ( row.child.isShown() ) {
            // This row is already open - close it
            $('div.slider', row.child()).slideUp( function () {
                row.child.hide();
                tr.removeClass('shown');
            } );
        }
        else {
            // Open this row
            row.child( createExpandedRowMarkup(row.data()), 'no-padding childRow' ).show();
            tr.addClass('shown');

            $('div.slider', row.child()).slideDown();
        }
    } );
    $(document.body).on('dblclick', '.childRow, .detailTable tbody tr td', function () {
        var table = $("#resultTable").DataTable();
        var tr = $(this).closest("tr.childRow").prev();
        var row = table.row( tr );
        var data = row.data();
        var toggleCell = tr.find(".toggle-details-row");
        // alert( 'You clicked on '+data['Target Gene Symbol']+'\'s row' );

        if ( row.child.isShown() ) {
            // This row is already open - close it
            $('div.slider', row.child()).slideUp( function () {
                row.child.hide();
                tr.removeClass('shown');
            } );
        }
        else {
            // Open this row
            row.child( createExpandedRowMarkup(row.data()), 'no-padding childRow' ).show();
            tr.addClass('shown');

            $('div.slider', row.child()).slideDown();
        }
    } );


    function displayAnalysis(analysis){
        analysis.displayInfo();
    }

    $("#download_result").click(function(){
        var dir_name = $("#result_file_select").val();
        console.log(dir_name);
        $("#download_result").text("Preparing .zip ...");
        $(location).attr('href', '/downloads/'+dir_name)
        $("#download_result").text("Download");
    });

    $("#display_settings_select").change( function() {
        old_display_type = display_type
        display_type = $(this).val();
        console.log(display_type);
        try {
            displayGeneBoxes(display_type);
        } catch(err) {
            alert("This analysis cannot be loaded for this file. "
            +"It's possible that this form of analysis was not performed on this dataset");
            console.log(err);
            $(this).val(old_display_type);
        }
    });

    $("#search").on("keyup", function(){
        //Provides search functionality for genes (filter on keyup)
        var input, filter, container, li, a, i;
        input = $("#search");
        filter = input.val().toUpperCase();
        container = $(".results");
        geneResults = container.children(".geneResult");
        for (i = 0; i < geneResults.length; i++) {
            a = $(geneResults[i]);
            if (a.attr('name').toUpperCase().indexOf(filter) > -1) {
                geneResults[i].style.display = "";
            } else {
                geneResults[i].style.display = "none";

            }
        }
    });
    $("div.results").on('click', '.geneResult', function(){
        //Expands gene box to display information
        if ($(this).hasClass("restrained")) {
            $(this).find(".details").removeClass("invisible");
            $(this).removeClass("restrained");
            console.log($(this));
            var $element = $( $(this) );
            setTimeout(function(){
                console.log($(this));
                var $window = $(window),
                    elementTop = $element.offset().top,
                    elementHeight = $element.height(),
                    viewportHeight = $window.height(),
                    scrollIt = elementTop - ((viewportHeight - elementHeight) / 2);
                $("HTML, BODY").animate({scrollTop: scrollIt}, 500);
            }, 200);
        } else {
            $(this).find(".details").addClass("invisible");
            $(this).addClass("restrained");
        }
        $('.geneResult').not(this).addClass("restrained");
        $('.geneResult').not(this).find(".details").addClass("invisible");
    });

    function scrollToCenter (my_element){
        console.log("Scroll?");
        var $window = $(window),
            $element = $(my_element),
            elementTop = $element.offset().top,
            elementHeight = $element.height(),
            viewportHeight = $window.height(),
            scrollIt = elementTop - ((viewportHeight - elementHeight) / 2);
            $("HTML, BODY").animate({scrollTop: scrollIt}, 500)
    };

});
