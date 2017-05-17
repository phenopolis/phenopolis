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
  PP.initLoginForm = function() {
    $('#login_modal').modal({ complete: function() { PP.showLoginForm(); } }); // Init Modals and ensure 
    PP.loginValidation();
    // Submit form upon clicking the login button in the login modal
    $('#login_btn').on('click', function() {
      $('#login_form').submit();
    });
    // submit form upon pressing enter in the login form
    $('#login_form input').keydown(function(e) {
      if (e.keyCode == 13) $('#login_form').submit();
    });
    // Open open Login Model On Clicking #enter_login_btn
    $('#enter_login_btn').on('click', function() {
      PP.showLoginForm();
    });
    PP.demoAutomaticLogin();
  };

  PP.loginValidation = function() {
    $('#login_form').validate({
      rules: {
        name: { 
          required: true
        },
        password: {
          required: true
        },
      },
      submitHandler: function(form) {
        $('#auth_modal').modal({ dismissible: false, endingTop: '20%' });
        $('#auth_modal').modal('open');

        $.ajax({
          type: 'POST',
          url: '/login',
          data: $('#login_form').serialize(),
          dataType: 'json',
          timeout: 120000,
          success: function(data) {
            window.location.href = '/search';
          },
          error: function(xhr, msg) {
            $('#auth_modal').modal('close');
            $("#username, #password").addClass("invalid");
            $("#username, #password").prop("aria-invalid", "true");
            $('#login_form_error_msg').show();
          }
        });
      }
    });
  };

  //
  PP.showLoginForm = function() {
    $('#login_content').show();
    $('#register_btn').show();
    $('#login_btn').show();
    $('#demo_login_btn').show();

    $('#register_content').hide();
    $('#enter_login_btn').hide();
    $('#submit_register_btn').hide();
  };

  // Automatically login when selecting the demo login button
  PP.demoAutomaticLogin = function() {
    $('#demo_login_btn, #parallax_demo_login_btn').on('click', function() {
      $('#name').val('demo');
      $('#password').val('demo123');
      $('#login_form').submit();
    });
  };

  PP.initRegistrationForm = function() {
    PP.registrationValidation();
    PP.showRegistrationFormOnClick();

    $('#submit_register_btn').on('click', function() {
      $('#register_form').submit();
    });
  };

  PP.registrationValidation = function() {
    $('#register_form').validate({
      rules: {
        name: { 
          required: true
        },
        email: {
          required: true
        },
        affliation: {
          required: true
        },
        "group[]": {
          required: true
        }
      },
      submitHandler: function(form) {
        $.ajax({
          type: 'POST',
          url: '/register',
          data: $('#register_form').serialize(),
          dataType: 'json',
          timeout: 120000,
          success: function(data) {
            $('#registering_msg_success').text(data.message);
            $('#registration_error').hide();
            $('#registration_success').show();
          },
          error: function(xhr, msg) {
            $('#registering_msg_error').text(xhr.responseJSON.message);
            $('#registration_success').hide();
            $('#registration_error').show();
          }
        });
      }
    });
  };

  //
  PP.showRegistrationFormOnClick = function() {
    $('#register_btn').on('click', function() {
      $('#login_content').hide();
      $('#register_btn').hide();

      $('#login_btn').hide();
      $('#demo_login_btn').hide();

      $('#register_content').show();
      $('#enter_login_btn').show();
      $('#submit_register_btn').show();
    });
  };

  // make the profile boxes on the about page the same height
  PP.equalizeAboutBoxes = function() {
    PP.equalizeDivHeight('.collection');
    $(window).resize(function() {
      $('.collection').css('height', 'auto');
      PP.equalizeDivHeight('.collection');
    });
  };

  PP.equalizeDivHeight = function(div) {
    if ($(div).length === 0) return;
    var maxHeight = 0;
    $(div).each(function() {
      if ($(this).height() > maxHeight) { maxHeight = $(this).height(); }
    });
    $(div).css('height', maxHeight);
  };
}());
// End of PP module
