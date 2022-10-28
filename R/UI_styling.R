# App styling css.
css <- tags$head(
    # Note the wrapping of the string in HTML()
    tags$style(HTML("
        .panel-body {
        padding: 5px;
        }
        .form-group {
        margin-bottom: 3px;
        padding-bottom: 2px !important;
        padding-top: 2px !important;
        font-size: 10px;
        line-height: 1.1;
        }
        .well {
        padding: 5px;
        margin-bottom: 10px;
        }
        .form-control, .selectize-input {
        padding-bottom: 2px !important;
        padding-top: 2px !important;
        font-size: 10px;
        height: 24px;
        min-height: 24px;
        line-height: 1.1;
        }
        .control-label {
        font-size: 10px;
        margin-bottom: 2px;
        }
        .panel-heading {
        padding: 5px 10px;
        }
        .selectize-control {
        margin-bottom: 0px;
        }
        .btn-file {
        padding: 1px 12px;
        }
        body {
        line-height: 1.1;
        }
        hr{
        margin-top: 3px;
        margin-bottom: 3px;
        }
        .panel-group {
        margin-bottom: 8px;
        }
        h3 {
        margin-top: 0px;
        }
        .nav li a.disabled {
        background-color: #aaa !important;
        color: #333 !important;
        cursor: not-allowed !important;
        border-color: #aaa !important;
        }
    "))
)