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
  PP.equalizeAboutBoxes = function() {
    PP.equalSizeDiv('.collection');
    $(window).resize(function() {
      $('.collection').css('height', 'auto');
      PP.equalSizeDiv('.collection');
    });
  };

  PP.equalSizeDiv = function(div) {
    if ($(div).length === 0) return;
    var maxHeight = 0;
    $(div).each(function() {
      if ($(this).height() > maxHeight) { maxHeight = $(this).height(); }
    });
    $(div).css('height', maxHeight);
  };

  PP.initLogin = function() {
    PP.openLoginModelOnClick('#enter_login_btn');

    // Submit form upon clicking the login button in the modal
    $('#login_btn').on('click', function() {
      PP.submitLogin();
    });

    // submit form upon pressing enter in the login form
    $('#login_form input').keydown(function(e) {
      if (e.keyCode == 13) PP.submitLogin();
    });

    // Automatically login when selecting the demo login button
    $('#demo_login_btn, #parallax_demo_login_btn').on('click', function() {
      $('#name').val('demo');
      $('#password').val('demo123');
      PP.submitLogin();
    });

    $('#register_btn').on('click', function() {
      PP.showRegistrationForm();
    });

    $('#submit_register_btn').on('click', function() {
      PP.submitRegisterForm();
    });
  };

  // TODO: Replace with JQuery.validate
  // show Login Model
  PP.submitLogin = function() {
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
      error: function(data, msg) {
        $('#auth_modal').modal('close');
        $("#username, #password").addClass("invalid");
        $("#username, #password").prop("aria-invalid", "true");
        $('#login_form_error_msg').show();
      }
    });
  };

  //
  PP.showRegistrationForm = function() {
    $('#login_content').hide();
    $('#register_btn').hide();
    $('#login_btn').hide();
    $('#demo_login_btn').hide();

    $('#register_content').show();
    $('#enter_login_btn').show();
    $('#submit_register_btn').show();
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

  //
  PP.openLoginModelOnClick = function(id) {
    $(id).on('click', function() {
      PP.showLoginForm();
    });
  };

  //
  PP.submitRegisterForm = function() {
    $.ajax({
      type: 'POST',
      url: '/register',
      data: $('#login_form').serialize(),
      dataType: 'json',
      timeout: 120000,
      success: function(data) {
        $('#registration_successful').show();
      },
      error: function(data, msg) {
        console.log(data);
        console.log(msg);
      }
    });
  };
}());
// End of PP module
