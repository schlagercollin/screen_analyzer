var global_dataframes = "Placeholder";
var display_type = "Placeholder";
var progress_counter = 0;

var analysis = "Placeholder";


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


function load_analysis(formData){
    // formData needs "analysis_name" key
    $.ajax({
        url: "/analysis/load",
        type: 'POST',
        data: formData,
        success: function(result) {
            $("button#load_result").html("<i class='fas fa-chart-line'></i> Load Analysis");
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
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "mageck_data pos_mageck"})
        } else if (val.startsWith("neg|")){
            formatted.push({"mData": val, "sTitle": val, "bVisible": false, "sClass": "mageck_data neg_mageck"})
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
        lengthMenu: [ 5, 10, 25, 50, 75, 100 ],
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
                                hide: ['.misc_data', '.count_data', '.pos_mageck', '.neg_mageck']
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
        load_analysis(formData);
    });

    //
    $(document.body).on('dblclick', '#resultTable tbody tr[role="row"]', function () {
        var table = $("#resultTable").DataTable();
        var tr = $(this);
        var row = table.row( tr );
        var data = row.data();
        var toggleCell = $(this).find(".toggle-details-row");
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

    function displayResult(output_dict) {

        //Parse output_dict and strings contained in there to create dataframes
        var nan = NaN; //needed for the eval of dataframe string
        for (var propt in output_dict) {
            obj_array = eval(output_dict[propt]);
            if (obj_array != false ) { //if there's data, set data in dataframe
                output_dict[propt] = jd.dfFromObjArray(obj_array);
            } else { //otherwise, keep field as false
                output_dict[propt] = false;
            }
        }

        global_dataframes = output_dict

        //Clear plotly plots that are present
        Plotly.purge('plotDiv_mageck');
        Plotly.purge('plotDiv_botCounts');
        Plotly.purge('plotDiv_unsortedCounts');
        Plotly.purge('plotDiv_topCounts');
        Plotly.purge('plotDiv_ratio');
        Plotly.purge('plotDiv_fischer');

        if (output_dict["mageck"]){ //if there's mageck data, extract x y and gene text data
            display_type = "Mageck";
            var x_mageck = output_dict["mageck"].c("pos|lfc").toArray();
            var genes_mageck = output_dict["mageck"].c("Target Gene Symbol").toArray();
            var y_mageck = output_dict["mageck"].c("-log(pos|p-value)").toArray();
        };

        if (output_dict["fischer"]){ // if there's fischer data, extract x y and gene text data
            display_type = "Fischer Exact"
            var x_fischer = output_dict["fischer"].c("LFC").toArray();
            var genes_fischer = output_dict["fischer"].c("Target Gene Symbol").toArray();
            var y_fischer = output_dict["fischer"].c("-log(FDR-Corrected P-Values)").toArray();
        };
        if (output_dict["ratio"]){ //if there's ratio test data, extract x y and gene text data
            display_type = "Ratio Test"
            var x_ratio = output_dict["ratio"].c("log_MA").toArray();
            var y_ratio = output_dict["ratio"].c("ZScore").toArray();
            var genes_ratio = output_dict["ratio"].c("Target Gene Symbol").toArray();
        };
        $("#display_settings_select").val(display_type);

        if (output_dict["top_sorted"]){
            var y_topCounts = output_dict["top_sorted"].c("Top Sorted Counts").toArray();
            var genes_topCounts = output_dict["top_sorted"].c("Target Gene Symbol").toArray();
            var x_topCounts = output_dict["top_sorted"].c("Target Gene Symbol").toArray();
        };
        if (output_dict["bot_sorted"]){
            var y_botCounts = output_dict["bot_sorted"].c("Bot Sorted Counts").toArray();
            var genes_botCounts = output_dict["bot_sorted"].c("Target Gene Symbol").toArray();
            var x_botCounts = output_dict["bot_sorted"].c("Target Gene Symbol").toArray();
        };
        if (output_dict["unsorted"]){
            var y_unsortedCounts = output_dict["unsorted"].c("Unsorted Counts").toArray();
            var genes_topCounts = output_dict["unsorted"].c("Target Gene Symbol").toArray();
            var x_unsortedCounts = output_dict["unsorted"].c("Target Gene Symbol").toArray();
        };

        // For a given display_type, display the geneboxes
        displayGeneBoxes(display_type)

        //Plotly stuff
        var trace_fischer = {
            x: x_fischer,
            y: y_fischer,
            text: genes_fischer,
            mode: 'markers',
            marker: {
                //color: colors,
                line: {width: 1}
            }
        }
        var trace_topCounts = {
            x: x_topCounts,
            y: y_topCounts,
            mode: 'markers',
            marker: {
                //color: colors,
                line: {width: 1}
            }
        }
        var trace_unsortedCounts = {
            x: x_unsortedCounts,
            y: y_unsortedCounts,
            mode: 'markers',
            marker: {
                //color: colors,
                line: {width: 1}
            }
        }
        var trace_mageck = {
            x: x_mageck,
            y: y_mageck,
            mode: 'markers',
            text: genes_mageck,
            marker: {
                //color: colors,
                line: {width: 1}
            }
        }
        var trace_botCounts = {
            x: x_botCounts,
            y: y_botCounts,
            text: genes_botCounts,
            mode: 'markers',
            marker: {
                //color: colors,
                line: {width: 1}
            }
        }
        var trace_ratio = {
            x: x_ratio,
            y: y_ratio,
            text: genes_ratio,
            mode: 'markers',
            marker: {
                //color: colors,
                line: {width: 1}
            }
        }
        var data_fischer = [ trace_fischer ];
        var data_topCounts = [ trace_topCounts ];
        var data_unsortedCounts = [ trace_unsortedCounts ];
        var data_mageck = [ trace_mageck ];
        var data_botCounts = [ trace_botCounts ];
        var data_ratio = [ trace_ratio ];

        var layout_fischer = {
            autosize: true,
            //width: 1300,
            //height: 600,
            title: "Fischer Test – log(P-Value) vs Log2(Fold Change))",
            xaxis: {
                autorange: true,
                title: "log2(FoldChange)"
            },
            yaxis: {
                //type: 'log',
                autorange: true,
                title: "-log(P-value)"
            },
            hovermode:'closest'
        }
        var layout_topCounts = {
            autosize: true,
            // width: 1300,
            // height: 600,
            title: "Top Sorted Frequencies",
            xaxis: {
                autorange: true,
                title: "Gene Target"
            },
            yaxis: {
                type: 'log',
                autorange: true,
                title: "Frequency"
            },
            hovermode:'closest'
        }
        var layout_unsortedCounts = {
            autosize: true,
            // width: 1300,
            // height: 600,
            title: "Unsorted Frequencies",
            xaxis: {
                autorange: true,
                title: "Gene Target"
            },
            yaxis: {
                type: 'log',
                autorange: true,
                title: "Frequency"
            },
            hovermode:'closest'
        }
        var layout_mageck = {
            autosize: true,
            // width: 1300,
            // height: 600,
            title: "Mageck – log(P-Value) vs log2(FoldChange)",
            xaxis: {
                autorange: true,
                title: "log2(FoldChange)"
            },
            yaxis: {
                //type: 'log',
                autorange: true,
                title: "-log(P-value)"
            },
            hovermode:'closest'
        }
        var layout_botCounts = {
            autosize: true,
            // width: 1300,
            // height: 600,
            title: "Bot Sorted Frequencies",
            xaxis: {
                autorange: true,
                title: "Gene Target"
            },
            yaxis: {
                //type: 'log',
                autorange: true,
                title: "Frequency"
            },
            hovermode:'closest'
        }
        var layout_ratio = {
            autosize: true,
            //width: 1300,
            //height: 600,
            title: "Regev Ratio Z Score",
            xaxis: {
                autorange: true,
                title: "log(Mean Abundance)"
            },
            yaxis: {
                //type: 'log',
                autorange: true,
                title: "Z Score"
            },
            hovermode:'closest'
        }

        $(".result_wrapper").fadeIn(100); //fade in result_wrapper to display results

        //If data for a given analyis type, create a plotly plot for it
        if (output_dict["fischer"]){
            Plotly.newPlot('plotDiv_fischer', data_fischer, layout_fischer); //create plot
        };
        if (output_dict["mageck"]){
            Plotly.newPlot('plotDiv_mageck', data_mageck, layout_mageck); //create plot
        };
        if (output_dict["ratio"]){
            Plotly.newPlot('plotDiv_ratio', data_ratio, layout_ratio); //create plot
        };
        if (output_dict["top_sorted"]){
            Plotly.newPlot('plotDiv_topCounts', data_topCounts, layout_topCounts); //create plot
        };
        if (output_dict["bot_sorted"]){
            Plotly.newPlot('plotDiv_botCounts', data_botCounts, layout_botCounts); //create plot
        };
        if (output_dict["unsorted"]){
            Plotly.newPlot('plotDiv_unsortedCounts', data_unsortedCounts, layout_unsortedCounts); //create plot
        };
    };

    function displayGeneBoxes(datatype){
        // Depending on datatype selection, print gene boxes and the information contained in them
        if (datatype == "Mageck") {
            console.log("Printing mageck boxes.");
            dataframe = global_dataframes["mageck"]
            genes_object = dataframe.toObjArray()
            var counter = 0;
            $("div.results").empty(); //clear old data stuff
            for (var gene in genes_object) {
                //Iterates through genes and displays information
                //HTML for each gene box
                gene_symbol = genes_object[gene]["Target Gene Symbol"]
                gene_id = genes_object[gene]["Target Gene ID"].toString()
                log2foldchange = parseFloat(genes_object[gene]["pos|lfc"]).toFixed(5);
                p_value = parseFloat(genes_object[gene]["-log(pos|p-value)"]).toFixed(5);
                unsorted_freq = parseInt(genes_object[gene]["Unsorted Counts"]).toFixed(5);
                top_sorted_freq = parseInt(genes_object[gene]["Top Sorted Counts"]).toFixed(5);
                bot_sorted_freq = parseInt(genes_object[gene]["Bot Sorted Counts"]).toFixed(5);
                summary = genes_object[gene]["Summary"]
                description = genes_object[gene]["Description"]
                //console.log(typeof(summary));
                if (isNaN(summary)){
                    //console.log("Replacing");
                    summary = "";
                }

                $("div.results").append(
                    "<div class='flex-item geneResult restrained' name="+genes_object[gene]["Target Gene Symbol"]+"><p><b>"+genes_object[gene]["Target Gene Symbol"]+"</b></p>" +
                    "<p title='log2(FoldChange)' style='font-size:0.5em'>"+log2foldchange+"</p>" +
                    "<p title='-log10(P-Value)' style='font-size:0.5em'>"+p_value+"</p>" +
                    "<p class='details invisible' style='font-size:0.5em'>"+description+"</p>" +
                    "<a class='no_style_link details invisible' style='font-size: 0.5em' href='https://www.ncbi.nlm.nih.gov/gene/"+gene_id+"'>See NCBI</a>"+
                    "<p class='details invisible' style='font-size:0.3em'>"+summary+"</p>" +
                    "</div>"
                );
                counter += 1;
                // if (counter >= 1000){
                //     break;
                // }
            };
        } else if (datatype=="Fischer Exact") {
            console.log("Printing fischer boxes");
            dataframe = global_dataframes["fischer"];
            genes_object = dataframe.toObjArray()
            var counter = 0;
            $("div.results").empty(); //clear old data stuff
            for (var gene in genes_object) {
                //Iterates through genes and displays information
                //HTML for each gene box
                gene_symbol = genes_object[gene]["Target Gene Symbol"]
                gene_id = genes_object[gene]["Target Gene ID"].toString()
                log2foldchange = parseFloat(genes_object[gene]["LFC"]).toFixed(5);
                p_value = parseFloat(genes_object[gene]["-log(FDR-Corrected P-Values)"]).toFixed(5);
                unsorted_freq = parseInt(genes_object[gene]["Unsorted Counts"]).toFixed(5);
                sorted_freq = parseInt(genes_object[gene]["Sorted Counts"]).toFixed(5);
                summary = genes_object[gene]["Summary"]
                description = genes_object[gene]["Description"]
                //console.log(typeof(summary));
                if (isNaN(summary)){
                    //console.log("Replacing");
                    summary = "";
                }

                $("div.results").append(
                    "<div class='flex-item geneResult restrained' name="+genes_object[gene]["Target Gene Symbol"]+"><p><b>"+genes_object[gene]["Target Gene Symbol"]+"</b></p>" +
                    "<p style='font-size:0.5em'>"+log2foldchange+"</p>" +
                    "<p style='font-size:0.5em'>"+p_value+"</p>" +
                    "<p class='details invisible' style='font-size:0.5em'>"+description+"</p>" +
                    "<a class='no_style_link details invisible' style='font-size: 0.5em' href='https://www.ncbi.nlm.nih.gov/gene/"+gene_id+"'>See NCBI</a>"+
                    "<p class='details invisible' style='font-size:0.3em'>"+summary+"</p>" +
                    "</div>"
                );
                counter += 1;
                // if (counter >= 1000){
                //     break;
                // }
            };
        } else if (datatype=="Ratio Test") {
            console.log("Printing ratio boxes");
            dataframe = global_dataframes["ratio"];
            console.log("toObjArray");
            genes_object = dataframe.toObjArray()
            var counter = 0;
            $("div.results").empty(); //clear old data stuff
            for (var gene in genes_object) {
                //Iterates through genes and displays information
                //HTML for each gene box
                gene_symbol = genes_object[gene]["Target Gene Symbol"]
                gene_id = genes_object[gene]["Target Gene ID"].toString()
                zscore = parseFloat(genes_object[gene]["ZScore"]).toFixed(5);
                log_MA = parseFloat(genes_object[gene]["log_MA"]).toFixed(5);
                top_sorted_freq = parseInt(genes_object[gene]["Top Sorted Counts"]).toFixed(5);
                bot_sorted_freq = parseInt(genes_object[gene]["Bot Sorted Counts"]).toFixed(5);
                summary = genes_object[gene]["Summary"]
                description = genes_object[gene]["Description"]
                //console.log(typeof(summary));
                if (isNaN(summary)){
                    //console.log("Replacing");
                    summary = "";
                }

                $("div.results").append(
                    "<div class='flex-item geneResult restrained' name="+genes_object[gene]["Target Gene Symbol"]+"><p><b>"+genes_object[gene]["Target Gene Symbol"]+"</b></p>" +
                    "<p title='Z Score' style='font-size:0.5em'>"+zscore+"</p>" +
                    "<p title='log_MA' style='font-size:0.5em'>"+log_MA+"</p>" +
                    "<p class='details invisible' style='font-size:0.5em'>"+description+"</p>" +
                    "<a class='no_style_link details invisible' style='font-size: 0.5em' target='_blank' href='https://www.ncbi.nlm.nih.gov/gene/"+gene_id+"'>See NCBI</a>"+
                    "<p class='details invisible' style='font-size:0.3em'>"+summary+"</p>" +
                    "</div>"
                );
                counter += 1;
                // console.log(counter);
                // if (counter >= 1000){
                //     break;
                // }
            };
        } else if (datatype=="Unsorted Counts") {
            console.log("Printing unsorted counts boxes");
            dataframe = global_dataframes["unsorted"];
            console.log("toObjArray");
            genes_object = dataframe.toObjArray()
            var counter = 0;
            $("div.results").empty(); //clear old data stuff
            for (var gene in genes_object) {
                //Iterates through genes and displays information
                //HTML for each gene box
                gene_symbol = genes_object[gene]["Target Gene Symbol"]
                gene_id = genes_object[gene]["Target Gene ID"].toString()
                unsorted_counts = parseInt(genes_object[gene]["Unsorted Counts"])
                top_sorted_freq = parseInt(genes_object[gene]["Top Sorted Counts"])
                bot_sorted_freq = parseInt(genes_object[gene]["Bot Sorted Counts"])
                summary = genes_object[gene]["Summary"]
                description = genes_object[gene]["Description"]
                //console.log(typeof(summary));
                if (isNaN(summary)){
                    //console.log("Replacing");
                    summary = "";
                }

                $("div.results").append(
                    "<div class='flex-item geneResult restrained' name="+genes_object[gene]["Target Gene Symbol"]+"><p><b>"+genes_object[gene]["Target Gene Symbol"]+"</b></p>" +
                    "<p title='Unsorted Counts' style='font-size:0.5em'>"+unsorted_counts+"</p>" +
                    "<p title='Top Sorted Counts' style='font-size:0.5em'>"+top_sorted_freq+"</p>" +
                    "<p title='Bot Sorted Counts' style='font-size:0.5em'>"+bot_sorted_freq+"</p>" +
                    "<p class='details invisible' style='font-size:0.5em'>"+description+"</p>" +
                    "<a class='no_style_link details invisible' style='font-size: 0.5em' target='_blank' href='https://www.ncbi.nlm.nih.gov/gene/"+gene_id+"'>See NCBI</a>"+
                    "<p class='details invisible' style='font-size:0.3em'>"+summary+"</p>" +
                    "</div>"
                );
                counter += 1;
                // console.log(counter);
                // if (counter >= 1000){
                //     break;
                // }
            };
        } else if (datatype=="Sorted or Top Sorted Counts") {
            console.log("Printing unsorted counts boxes");
            dataframe = global_dataframes["top_sorted"];
            console.log("toObjArray");
            genes_object = dataframe.toObjArray()
            var counter = 0;
            $("div.results").empty(); //clear old data stuff
            for (var gene in genes_object) {
                //Iterates through genes and displays information
                //HTML for each gene box
                gene_symbol = genes_object[gene]["Target Gene Symbol"]
                gene_id = genes_object[gene]["Target Gene ID"].toString()
                unsorted_counts = parseInt(genes_object[gene]["Unsorted Counts"])
                top_sorted_freq = parseInt(genes_object[gene]["Top Sorted Counts"])
                bot_sorted_freq = parseInt(genes_object[gene]["Bot Sorted Counts"])
                summary = genes_object[gene]["Summary"]
                description = genes_object[gene]["Description"]
                //console.log(typeof(summary));
                if (isNaN(summary)){
                    //console.log("Replacing");
                    summary = "";
                }

                $("div.results").append(
                    "<div class='flex-item geneResult restrained' name="+genes_object[gene]["Target Gene Symbol"]+"><p><b>"+genes_object[gene]["Target Gene Symbol"]+"</b></p>" +
                    "<p title='Top Sorted Counts' style='font-size:0.5em'>"+top_sorted_freq+"</p>" +
                    "<p title='Bot Sorted Counts' style='font-size:0.5em'>"+bot_sorted_freq+"</p>" +
                    "<p title='Unsorted Counts' style='font-size:0.5em'>"+unsorted_counts+"</p>" +

                    "<p class='details invisible' style='font-size:0.5em'>"+description+"</p>" +
                    "<a class='no_style_link details invisible' style='font-size: 0.5em' target='_blank' href='https://www.ncbi.nlm.nih.gov/gene/"+gene_id+"'>See NCBI</a>"+
                    "<p class='details invisible' style='font-size:0.3em'>"+summary+"</p>" +
                    "</div>"
                );
                counter += 1;
                // console.log(counter);
                // if (counter >= 1000){
                //     break;
                // }
            };
        } else if (datatype=="Bottom Sorted Counts") {
            console.log("Printing unsorted counts boxes");
            dataframe = global_dataframes["bot_sorted"];
            console.log("toObjArray");
            genes_object = dataframe.toObjArray()
            var counter = 0;
            $("div.results").empty(); //clear old data stuff
            for (var gene in genes_object) {
                //Iterates through genes and displays information
                //HTML for each gene box
                gene_symbol = genes_object[gene]["Target Gene Symbol"]
                gene_id = genes_object[gene]["Target Gene ID"].toString()
                unsorted_counts = parseInt(genes_object[gene]["Unsorted Counts"])
                top_sorted_freq = parseInt(genes_object[gene]["Top Sorted Counts"])
                bot_sorted_freq = parseInt(genes_object[gene]["Bot Sorted Counts"])
                summary = genes_object[gene]["Summary"]
                description = genes_object[gene]["Description"]
                //console.log(typeof(summary));
                if (isNaN(summary)){
                    //console.log("Replacing");
                    summary = "";
                }

                $("div.results").append(
                    "<div class='flex-item geneResult restrained' name="+genes_object[gene]["Target Gene Symbol"]+"><p><b>"+genes_object[gene]["Target Gene Symbol"]+"</b></p>" +
                    "<p title='Bot Sorted Counts' style='font-size:0.5em'>"+bot_sorted_freq+"</p>" +
                    "<p title='Top Sorted Counts' style='font-size:0.5em'>"+top_sorted_freq+"</p>" +
                    "<p title='Unsorted Counts' style='font-size:0.5em'>"+unsorted_counts+"</p>" +

                    "<p class='details invisible' style='font-size:0.5em'>"+description+"</p>" +
                    "<a class='no_style_link details invisible' style='font-size: 0.5em' target='_blank' href='https://www.ncbi.nlm.nih.gov/gene/"+gene_id+"'>See NCBI</a>"+
                    "<p class='details invisible' style='font-size:0.3em'>"+summary+"</p>" +
                    "</div>"
                );
                counter += 1;
                // console.log(counter);
                // if (counter >= 1000){
                //     break;
                // }
            };
        } else {
            console.log("No data type selected.");
        };

    }

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
