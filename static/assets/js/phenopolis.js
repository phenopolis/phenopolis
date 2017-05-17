/*
    PP - Phenopolis's JavaScript module

    Define a global PP (acronym for Phenopolis) object containing all
    PP associated methods:
*/

// define global PP object
var PP;
if (!PP) {
  PP = {};
}

// PP module
(function() {
  // Adds the side navbar dropdown
  PP.addUserDropDown = function() {
    $('.user .dropdown-button').dropdown({
      inDuration: 300,
      outDuration: 225,
      hover: true,
      belowOrigin: true,
      alignment: 'right'
    });
  };

  /*
   * AUTOCOMPLETE FUNCTIONS 
   */

  // Sets up Autocomplete & Allows submission upon pressing enter
  PP.setUpSearchField = function() {
    var inputAutoComplete = PP.initialiseBloodHound('/autocomplete');
    PP.initialiseTypeahead('.searching .typeahead', inputAutoComplete);
    PP.fixHomeSearchBug();
    PP.addNavbarSearchAnimation();
    PP.bindTypeaheadSelect();
    PP.submitSearchOnEnter();
  };

  // Used in autocomplete - connects to the server and provides autocomplete options
  PP.initialiseBloodHound = function(autocompleteUrl) {
    return new Bloodhound({
      datumTokenizer: Bloodhound.tokenizers.obj.whitespace('value'),
      queryTokenizer: Bloodhound.tokenizers.whitespace,
      remote: {
        url: autocompleteUrl + '/%QUERY',
        wildcard: '%QUERY'
      }
    });
  };

  //
  PP.initialiseTypeahead = function(input, inputAutoComplete) {
    $(input).typeahead({
      autoselect: true,
      classNames: { menu: "autocomplete-content dropdown-content" },
    }, {
      name: 'my-dataset',
      source: inputAutoComplete,
      limit: 19
    });

  };



  // Typeahead moves the input into a div. This causes issues with other JS/CSS
  // because the label and icon are no longer sister to the input
  // So simply move these into the same div as the input...
  PP.fixHomeSearchBug = function() {
    $('#search_home.typeahead').after($('#home_search_group label'));
    $('#search_home.typeahead').before($('#home_search_group .material-icons'));
  };

  //
  PP.addNavbarSearchAnimation = function() {
    $("#navbar_search_form").on('focus', '#search_navbar', function() {
      $('#search_navbar_label .material-icons').css('color', '#444');
      $('#search_navbar_cross').css('color', '#444');
      $('#navbar_search_form').css('background-color', '#fff');
    });
    $("#navbar_search_form").on('blur', '#search_navbar', function() {
      $('#search_navbar_label .material-icons').css('color', '#fff');
      $('#search_navbar_cross').css('color', 'transparent');
      $('#navbar_search_form').css('background-color', 'rgba(79, 166, 214, 0.5)');
    });
  };

  // redirect page upon selecting an option
  PP.bindTypeaheadSelect = function() {
    $('.typeahead').bind('typeahead:select', function(obj, datum) {
      window.location.href = '/awesome?query=' + datum;
    });
  };

  // Allows submission on enter with the autocomplete search bar
  PP.submitSearchOnEnter = function() {
    $('.typeahead').keypress(function(e) {
      if (e.which == 13) {
        window.location.href = '/awesome?query=' + $(this).val();
      }
    });
  };

   PP.setUpValidatorDefaults = function() {
    $.validator.addMethod("secure_password", function(value, element) {
      return /^[A-Za-z0-9\d=!\-@._*]*$/.test(value) && // consists of only these 
        /[a-z]/.test(value) && // has a lowercase letter
        /[A-Z]/.test(value) && // has an uppercase letter
        /\d/.test(value); // has a digit
    }, 'Password must contain an upper case letter, a lower case letter and a digit.');

    $.validator.setDefaults({
      errorClass: 'invalid',
      validClass: "valid",
      errorPlacement: function (error, element) {
          $(element)
            .closest("form")
            .find("label[for='" + element.attr("id") + "']")
            .attr('data-error', error.text());
      },
    });
  };

  PP.initChangePasswordForm = function() {
    PP.changePasswordValidation();

    $('#submit_change_password_btn').on('click', function() {
      $('#change_password_form').submit();
    });
  };

  PP.changePasswordValidation = function() {
    $('#change_password_form').validate({
      rules: {
        current_password: { 
          required: true,
        },
        new_password_1: {
          required: true,
          secure_password: true,
          minlength: 6
        },
        new_password_2: {
          required: true,
          secure_password: true,
          minlength: 6,
          equalTo: '#new_password_1'
        },
      },
      messages: {
        new_password_2: {
          equalTo: 'Both new passwords must match.'
        }
      },
      submitHandler: function(form) {
        $('#auth_modal').modal({ dismissible: false, endingTop: '20%' });
        $('#auth_modal').modal('open');
        $('#change_password_form_error_msg').hide();
        $('#change_password_successful').hide();

        $.ajax({
            type: 'POST',
            url: '/change_password',
            data: $('#change_password_form').serialize(),
            dataType: 'json',
            timeout: 120000,
            success: function (data) {
                $('#auth_modal').modal('close');
                $('#change_password_successful').show();
                $("#change_password_successful").text(data.success);
            },
            error: function (data, msg) {
                $('#auth_modal').modal('close');
                $("#username, #password, #new_password_1, #new_password_2").addClass("invalid");
                $("#username, #password, #new_password_1, #new_password_2").prop("aria-invalid", "true");
                $('#change_password_form_error_msg').show();
                $("#change_password_form_error_msg").text(data.responseJSON.error);
            }
        });
      }
    });
  };

  /*
   *  Generic  Functions
   */

  PP.create_url = function(base_url, text, url_end, a_class) {
    if (url_end === undefined) {
      return '<a href=' + base_url + '/' + text + ' target="_blank">' + text + '</a>';
    } else if (a_class === undefined) {
      return '<a href=' + base_url + '/' + url_end + ' target="_blank">' + text + '</a>';
    } else {
      return '<a class="' + a_class + '" href=' + base_url + '/' + url_end + ' target="_blank">' + text + '</a>';
    }
  };

  PP.initTableSorter = function(table) {
    $(table).tablesorter({
      theme: 'materialize',
      widthFixed: true, // hidden filter input/selects will resize the columns, so try to minimize the change
      widgets: ["filter"], // initialize zebra striping and filter widgets
      widgetOptions: {
        // Use the $.tablesorter.storage utility to save the most recent filters
        filter_saveFilters: true,
        // jQuery selector string of an element used to reset the filters
        filter_reset: '.reset',
        // extra css class name (string or array) added to the filter element (input or select)
        // select needs a "browser-default" class or it gets hidden
        filter_placeholder: { search: 'Filter' },
        filter_columnFilters: true,
        filter_matchType : { 'input': 'match', 'select': 'match' }
      }
    });
    $(table).tablesorterPager({
      container: $(table+" .ts-pager"),
      storageKey: table.replace('#','')+'pager',
      // // remove rows from the table to speed up the sort of large tables.
      // // setting this to false, only hides the non-visible rows; needed if you plan to add/remove rows with the pager enabled.
      removeRows: true,

      output: '{startRow}-{endRow} of {filteredRows}  ({totalRows})'
    });
    $('td:has(.chip)').css('text-align', 'left');
    $(table).bind('pagerComplete', function() {
      $('td:has(.chip)').css('text-align', 'left');
    });
  };
}());
// End of PP module
