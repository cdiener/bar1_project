#include <gtk/gtk.h>


gboolean
on_segmentDialog_delete_event          (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_cbAlignCells_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbAlignBoundaries_toggled           (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbNucleus_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbRecomb_toggled                    (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_eCell_editing_done                  (GtkCellEditable *celleditable,
                                        gpointer         user_data);

void
on_eNucleus_editing_done               (GtkCellEditable *celleditable,
                                        gpointer         user_data);

void
on_eNew_editing_done                   (GtkCellEditable *celleditable,
                                        gpointer         user_data);

void
on_cboxFlagCell_changed                (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cboxFlagNucleus_changed             (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cboxFlagNew_changed                 (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cancelbuttonSeg_clicked             (GtkButton       *button,
                                        gpointer         user_data);

void
on_okbuttonSeg_clicked                 (GtkButton       *button,
                                        gpointer         user_data);

void
on_applybuttonSeg_clicked              (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_app1_delete_event                   (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_new1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_open1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_save1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_save_as1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_quit1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_image_lists1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_output_config1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_segmentation1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_default_segm_val1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_run_test1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_analyze1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_single_image1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_single_pos2_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_batch1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_browse1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_about1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_about2_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_tbInput_clicked                     (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_tbSegmentation_clicked              (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_tbTest_clicked                      (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

gboolean
on_inputDialog_delete_event            (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_fcbBrightField_selection_changed    (GtkFileChooser  *filechooser,
                                        gpointer         user_data);

void
on_fcbFluo_selection_changed           (GtkFileChooser  *filechooser,
                                        gpointer         user_data);

void
on_fcbThird_selection_changed          (GtkFileChooser  *filechooser,
                                        gpointer         user_data);

void
on_fcbFlat_selection_changed           (GtkFileChooser  *filechooser,
                                        gpointer         user_data);

void
on_fcbDark_selection_changed           (GtkFileChooser  *filechooser,
                                        gpointer         user_data);

void
on_cbForcePath_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_fcbMapping_selection_changed        (GtkFileChooser  *filechooser,
                                        gpointer         user_data);

void
on_cboxFretBright_changed              (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cbFret_toggled                      (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cboxFretNucleus_changed             (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cboxImageType_changed               (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cboxMapping_changed                 (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cboxThird_changed                   (GtkComboBox     *combobox,
                                        gpointer         user_data);

void
on_cbBFasFL_toggled                    (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_rbNoAlignment_group_changed         (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_rbFirstAlignment_toggled            (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_rbBRAlignment_toggled               (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbIndividual_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cancelbuttonInput_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_okbuttonInput_clicked               (GtkButton       *button,
                                        gpointer         user_data);

void
on_applybuttonInput_clicked            (GtkButton       *button,
                                        gpointer         user_data);

void
on_rbAlignCells_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_rbNoAlignCells_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_rbAlignCellsBoundaries_toggled      (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbT0_toggled                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbNewRecomb_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbNucRecomb_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbFLRecomb_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_bcancelImageSetup_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_bApplyImageSetup_clicked            (GtkButton       *button,
                                        gpointer         user_data);

void
on_bOkImageSetup_clicked               (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_ImageSetupDialog_delete_event       (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_tbImageSetup_clicked                (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_output1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_cbPosToken_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbTimeToken_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbSeparator_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_entryFL_changed                     (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_entrySeparator_changed              (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_entryPosToken_changed               (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_entryTimeToken_changed              (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_inputDialog_show                    (GtkWidget       *widget,
                                        gpointer         user_data);

void
on_cbGradient_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbBackground_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_cbNucleus_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_tbZoomIn_clicked                    (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_tbZoomOut_clicked                   (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_tbZoom1_clicked                     (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_tbZoomFit_clicked                   (GtkToolButton   *toolbutton,
                                        gpointer         user_data);

void
on_treeview1_cursor_changed            (GtkTreeView     *treeview,
                                        gpointer         user_data);

gboolean
on_treeview1_popup_menu                (GtkWidget       *widget,
                                        gpointer         user_data);

gboolean
on_treeview1_button_press_event        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data);

gboolean
on_treeview1_button_press_event        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data);

gboolean
on_treeview1_popup_menu                (GtkWidget       *widget,
                                        gpointer         user_data);

void
on_cbTreeViewOutFilter_toggled         (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_bCancel_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_bOk_clicked                         (GtkButton       *button,
                                        gpointer         user_data);

void
on_gsetup_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

gboolean
on_setupDialog_delete_event            (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_bTest_clicked                       (GtkButton       *button,
                                        gpointer         user_data);

void
on_okbutton2_clicked                   (GtkButton       *button,
                                        gpointer         user_data);

gboolean
on_calculatingDialog_delete_event      (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

gboolean
on_calculatingDialog_delete_event      (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_bRun_clicked                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_loadimages_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_imagesetup_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_segmentation_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

gboolean
on_aboutdialog1_delete_event           (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_aboutdialog1_response               (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data);
