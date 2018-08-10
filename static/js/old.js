// Rudimentary side by side expansion

function createExpandedRowMarkup ( d ) {
    // `d` is the original data object for the row
    left_markup = "";
    right_markup = "";
    var add_to_left = true;
    for ( const [ key, value ] of Object.entries(d)) {
        markup = `
            <tr>
                <td> ${ key } </td>
                <td> ${ value } </td>
            </tr>
        `
        if (add_to_left) {
            left_markup+=markup;
            console.log("left");
        } else {
            right_markup+=markup;
            console.log("right");
        };
        add_to_left = !add_to_left
    };

    full_markup = `
        <table cellpadding="5" cellspacing="0" border="0" style="width: 40%; padding-left:50px; display: block; float: left;">
            ${ left_markup }
        </table>
        <table cellpadding="5" cellspacing="0" border="0" style="width: 40%; padding-left:50px; display: block; float: left;">
            ${ right_markup }
        </table>
    `
    return full_markup;
}
