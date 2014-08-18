#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <gtk/gtk.h>
#include <glib.h>
#include <glib/gstdio.h>

#include "callbacks.h"
#include "interface.h"
#include "support.h"

#include "main.h"

extern GtkWidget *app1;
extern GtkWidget *inputDialog;
extern GtkWidget *segmentDialog;
extern GtkWidget *imageSetupDialog;
extern GtkWidget *setupDialog;
extern GtkAnimView *view;
extern GtkTreeStore *treestore;
extern GtkWidget    *aboutDialog;

extern PARAMS p;

enum
{
  STRING_COL = 0,
  TIMETAG_COL,
  NUM_COLS
} ;

gboolean
on_segmentDialog_delete_event          (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(widget)));
  return TRUE;
}


void
on_cbAlignCells_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_cbAlignBoundaries_toggled           (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}



void
on_cbRecomb_toggled                    (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_eCell_editing_done                  (GtkCellEditable *celleditable,
                                        gpointer         user_data)
{

}


void
on_eNucleus_editing_done               (GtkCellEditable *celleditable,
                                        gpointer         user_data)
{

}


void
on_eNew_editing_done                   (GtkCellEditable *celleditable,
                                        gpointer         user_data)
{

}


void
on_cboxFlagCell_changed                (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

}


void
on_cboxFlagNucleus_changed             (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

}


void
on_cboxFlagNew_changed                 (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

}


void
on_cancelbuttonSeg_clicked             (GtkButton       *button,
                                        gpointer         user_data)
{
  Param2SegmentationDialog();
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


void
on_okbuttonSeg_clicked                 (GtkButton       *button,
                                        gpointer         user_data)
{
 if(!CheckSegmentationForm()) return;
 SegmentationDialog2Param();
 gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


void
on_applybuttonSeg_clicked              (GtkButton       *button,
                                        gpointer         user_data)
{
 CheckSegmentationForm();
 SegmentationDialog2Param();
}


gboolean
on_app1_delete_event                   (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  //Cleaning TMP stuff before leaving
  CleanTMPdir();
  gtk_exit(0); 
  return FALSE;
}


void
on_new1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_open1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_save1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_save_as1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_quit1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
  gtk_main_quit();
}


void
on_image_lists1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_output_config1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_segmentation1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_default_segm_val1_activate          (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_run_test1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_analyze1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_single_image1_activate              (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_single_pos2_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_batch1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_browse1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_about1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_about2_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
    gtk_widget_show(aboutDialog);
    return;
}


void
on_tbInput_clicked                     (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
    gtk_widget_show(inputDialog);
    return;
}


void
on_tbSegmentation_clicked              (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
    gtk_widget_show(segmentDialog);
    return;
}


void
on_tbTest_clicked                      (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{

}


gboolean
on_inputDialog_delete_event            (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(widget)));
  return TRUE;
}


void
on_fcbBrightField_selection_changed    (GtkFileChooser  *filechooser,
                                        gpointer         user_data)
{

}


void
on_fcbFluo_selection_changed           (GtkFileChooser  *filechooser,
                                        gpointer         user_data)
{

}


void
on_fcbThird_selection_changed          (GtkFileChooser  *filechooser,
                                        gpointer         user_data)
{

}


void
on_fcbFlat_selection_changed           (GtkFileChooser  *filechooser,
                                        gpointer         user_data)
{

}


void
on_fcbDark_selection_changed           (GtkFileChooser  *filechooser,
                                        gpointer         user_data)
{

}


void
on_cbForcePath_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 
 if(!gtk_toggle_button_get_active(togglebutton)){
  w = lookup_widget(inputDialog,"fcbFluo");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);

  w = lookup_widget(inputDialog,"fcbFlat");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);

  w = lookup_widget(inputDialog,"fcbDark");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);

 }else{
  w = lookup_widget(inputDialog,"fcbFluo");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);

  w = lookup_widget(inputDialog,"fcbFlat");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);

  w = lookup_widget(inputDialog,"fcbDark");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
 }
}


void
on_fcbMapping_selection_changed        (GtkFileChooser  *filechooser,
                                        gpointer         user_data)
{

}


void
on_cbFret_toggled                      (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
   
   if(gtk_toggle_button_get_active(togglebutton)){
    w = lookup_widget(imageSetupDialog,"cboxFretNucleus");
    //gtk_widget_set(GTK_WIDGET(w),"visible",TRUE,NULL);
    g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
  }else{
     w = lookup_widget(imageSetupDialog,"cboxFretNucleus");
     g_object_set(G_OBJECT(w),"visible",FALSE,NULL);
   }
}


void
on_cboxFretNucleus_changed             (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

}


void
on_cboxImageType_changed               (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

}


void
on_cboxMapping_changed                 (GtkComboBox     *combobox,
                                        gpointer         user_data)
{
 GtkWidget *w;
 
 w = lookup_widget(inputDialog,"fcbMapping");
 if(GetComboBoxItem(GTK_COMBO_BOX(combobox),NULL)==BF_FL_MAPPING_LIST){
  gtk_widget_show(GTK_WIDGET(w));
 }else{
  gtk_widget_hide(GTK_WIDGET(w));
 }
}


void
on_cboxThird_changed                   (GtkComboBox     *combobox,
                                        gpointer         user_data)
{

}


void
on_cbBFasFL_toggled                    (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_rbNoAlignment_group_changed         (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_rbFirstAlignment_toggled            (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_rbBRAlignment_toggled               (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_cbIndividual_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_cancelbuttonInput_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
  Param2InputDialog();
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


void
on_okbuttonInput_clicked               (GtkButton       *button,
                                        gpointer         user_data)
{ 

  if(!CheckInputForm()) return;
  InputDialog2Param();
  GtkWidget *w;

  if(!LoadImageNames()) return;
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));

  w = lookup_widget(segmentDialog,"hbTest");
  g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
//  w = lookup_widget(segmentDialog,"frameAlignement");
//  g_object_set(G_OBJECT(w),"visible",FALSE,NULL);  
}

/*
void
on_applybuttonInput_clicked            (GtkButton       *button,
                                        gpointer         user_data)
{
  CheckInputForm();
  InputDialog2Param();
}
*/

void
on_rbAlignCells_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_rbNoAlignCells_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_rbAlignCellsBoundaries_toggled      (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_cbT0_toggled                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 
 if(gtk_toggle_button_get_active(togglebutton)){
  w = lookup_widget(segmentDialog,"sbMaxSplit0");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
  w = lookup_widget(segmentDialog,"sbMaxDist0");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
 }else{
  w = lookup_widget(segmentDialog,"sbMaxSplit0");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
  w = lookup_widget(segmentDialog,"sbMaxDist0");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
 }
}


void
on_cbNewRecomb_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_cbNucRecomb_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 
 if(gtk_toggle_button_get_active(togglebutton)){
  w = lookup_widget(segmentDialog,"cboxNucRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
  w = lookup_widget(segmentDialog,"entryNucRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
 }else{
  w = lookup_widget(segmentDialog,"cboxNucRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
  w = lookup_widget(segmentDialog,"entryNucRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
 }
}


void
on_cbFLRecomb_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 
 if(gtk_toggle_button_get_active(togglebutton)){
  w = lookup_widget(segmentDialog,"cboxFLRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
  w = lookup_widget(segmentDialog,"entryFLRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),TRUE);
 }else{
  w = lookup_widget(segmentDialog,"cboxFLRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
  w = lookup_widget(segmentDialog,"entryFLRecomb");
  gtk_widget_set_sensitive(GTK_WIDGET(w),FALSE);
 }
}


void
on_bcancelImageSetup_clicked           (GtkButton       *button,
                                        gpointer         user_data)
{
  Param2ImageSetupDialog();
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


void
on_bApplyImageSetup_clicked            (GtkButton       *button,
                                        gpointer         user_data)
{
 CheckImageSetupForm();
 ImageSetupDialog2Param();
}


void
on_bOkImageSetup_clicked               (GtkButton       *button,
                                        gpointer         user_data)
{
 CheckImageSetupForm();
 ImageSetupDialog2Param();
 gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


gboolean
on_ImageSetupDialog_delete_event       (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(widget)));
  return TRUE;
}


void
on_tbImageSetup_clicked                (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  w = lookup_widget(segmentDialog,"frameAlignement");
  g_object_set(G_OBJECT(w),"visible",TRUE,NULL);
  gtk_widget_show(imageSetupDialog);
  return;
}


void
on_output1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{

}


void
on_cbPosToken_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  gboolean flag;
  
  w = lookup_widget(inputDialog,"entryPosToken");
  flag = gtk_toggle_button_get_active(togglebutton);
  gtk_widget_set_sensitive(GTK_WIDGET(w),flag);
  //p.bpos = flag;
  UpdateEntryFileName();
}


void
on_cbTimeToken_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  gboolean flag;
  
  w = lookup_widget(inputDialog,"entryTimeToken");
  flag = gtk_toggle_button_get_active(togglebutton);
  gtk_widget_set_sensitive(GTK_WIDGET(w),flag);
  //p.btime = flag;
  UpdateEntryFileName();
}


void
on_cbSeparator_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
  GtkWidget *w;
  gboolean flag;
  
  w = lookup_widget(inputDialog,"entrySeparator");
  flag = gtk_toggle_button_get_active(togglebutton);
  gtk_widget_set_sensitive(GTK_WIDGET(w),flag);
  //p.bsep = flag;
  UpdateEntryFileName();
}

void
on_entryFL_changed                     (GtkEditable     *editable,
                                        gpointer         user_data)
{
 UpdateEntryFileName();
}


void
on_entrySeparator_changed              (GtkEditable     *editable,
                                        gpointer         user_data)
{
 UpdateEntryFileName();
}


void
on_entryPosToken_changed               (GtkEditable     *editable,
                                        gpointer         user_data)
{
 UpdateEntryFileName();
}


void
on_entryTimeToken_changed              (GtkEditable     *editable,
                                        gpointer         user_data)
{
 UpdateEntryFileName();
}


void
on_inputDialog_show                    (GtkWidget       *widget,
                                        gpointer         user_data)
{
 UpdateEntryFileName();
}


void
on_cbGradient_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 gboolean flag = gtk_toggle_button_get_active(togglebutton);
  
 w    = lookup_widget(inputDialog,"entryFlat");
 gtk_widget_set_sensitive(GTK_WIDGET(w),flag);
}


void
on_cbBackground_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 gboolean flag = gtk_toggle_button_get_active(togglebutton);
  
 w    = lookup_widget(inputDialog,"entryDark");
 gtk_widget_set_sensitive(GTK_WIDGET(w),flag);
}


void
on_cbNucleus_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{
 GtkWidget *w;
 gboolean flag = gtk_toggle_button_get_active(togglebutton);
  
  w    = lookup_widget(imageSetupDialog,"entryNucleus");
  gtk_widget_set_sensitive(GTK_WIDGET(w),flag);
}


void
on_tbZoomIn_clicked                    (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
	gtk_image_view_zoom_in (GTK_IMAGE_VIEW (view));
}


void
on_tbZoomOut_clicked                   (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
	gtk_image_view_zoom_out (GTK_IMAGE_VIEW (view));
}


void
on_tbZoom1_clicked                     (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
	gtk_image_view_set_zoom (GTK_IMAGE_VIEW (view), 1.0);
}


void
on_tbZoomFit_clicked                   (GtkToolButton   *toolbutton,
                                        gpointer         user_data)
{
	gtk_image_view_set_fitting (GTK_IMAGE_VIEW (view), TRUE);
}


void
on_treeview1_cursor_changed            (GtkTreeView     *treeview,
                                        gpointer         user_data)
{
  GtkTreeSelection *selection;
  GtkTreeIter       iter;

  /* This will only work in single or browse selection mode! */
  selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));
  if (gtk_tree_selection_get_selected(selection, NULL, &iter))
  {
    gchar *name;

    gtk_tree_model_get (GTK_TREE_MODEL(treestore), &iter, STRING_COL, &name, -1);
    if(g_strrstr(name,".tif") || g_strrstr(name,".TIF"))
     LoadImages(name);
    g_free(name);
  }
  else
  {
    g_print ("no row selected.\n");
  }

}


gboolean
on_treeview1_button_press_event        (GtkWidget       *widget,
                                        GdkEventButton  *event,
                                        gpointer         user_data)
{

  return FALSE;
}


gboolean
on_treeview1_popup_menu                (GtkWidget       *widget,
                                        gpointer         user_data)
{
  return FALSE;
}

/////////////////  

gboolean
  view_onButtonPressed (GtkWidget *treeview, GdkEventButton *event, gpointer userdata)
  {
    /* single click with the right mouse button? */
    if (event->type == GDK_BUTTON_PRESS  &&  event->button == 3)
    {
      g_print ("Single right click on the tree view.\n");

      /* optional: select row if no row is selected or only
       *  one other row is selected (will only do something
       *  if you set a tree selection mode as described later
       *  in the tutorial) */
      if (1)
      {
        GtkTreeSelection *selection;

        selection = gtk_tree_view_get_selection(GTK_TREE_VIEW(treeview));

        /* Note: gtk_tree_selection_count_selected_rows() does not
         *   exist in gtk+-2.0, only in gtk+ >= v2.2 ! */
        if (gtk_tree_selection_count_selected_rows(selection)  <= 1)
        {
           GtkTreePath *path;

           /* Get tree path for row that was clicked */
           if (gtk_tree_view_get_path_at_pos(GTK_TREE_VIEW(treeview),
                                             (gint) event->x, 
                                             (gint) event->y,
                                             &path, NULL, NULL, NULL))
           {
             gtk_tree_selection_unselect_all(selection);
             gtk_tree_selection_select_path(selection, path);
             gtk_tree_path_free(path);
           }
        }
      } /* end of optional bit */

      view_popup_menu(treeview, event, userdata);

      return TRUE; /* we handled this */
    }

    return FALSE; /* we did not handle this */
  }



/////////////////////////


void
on_cbTreeViewOutFilter_toggled         (GtkToggleButton *togglebutton,
                                        gpointer         user_data)
{

}


void
on_bCancel_clicked                     (GtkButton       *button,
                                        gpointer         user_data)
{
 Param2SetupDialog();
 gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


void
on_bOk_clicked                         (GtkButton       *button,
                                        gpointer         user_data)
{
 if(!CheckExecFile()) return;
 SetupDialog2Param();
 gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


void
on_gsetup_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
 gtk_widget_show(GTK_WIDGET(setupDialog));
}


gboolean
on_setupDialog_delete_event            (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(widget)));
  return TRUE;
}


void
on_bTest_clicked                       (GtkButton       *button,
                                        gpointer         user_data)
{
 if(!CheckSegmentationForm()) return;
 SegmentationDialog2Param();
 TestCellID();
}


void
on_okbutton2_clicked                   (GtkButton       *button,
                                        gpointer         user_data)
{
 gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(button)));
}


gboolean
on_calculatingDialog_delete_event      (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  gtk_widget_hide(gtk_widget_get_toplevel(GTK_WIDGET(widget)));
  return TRUE;
}


void
on_bRun_clicked                        (GtkButton       *button,
                                        gpointer         user_data)
{
 if(!CheckSegmentationForm()) return;
 SegmentationDialog2Param();
 RunCellID();
}


void
on_loadimages_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
    gtk_widget_show(inputDialog);
    return;
}


void
on_imagesetup_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
    gtk_widget_show(imageSetupDialog);
    return;
}


void
on_segmentation_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data)
{
    gtk_widget_show(segmentDialog);
    return;
}


gboolean
on_aboutdialog1_delete_event           (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data)
{
  gtk_widget_hide(aboutDialog);
  return TRUE;
}


void
on_aboutdialog1_response               (GtkDialog       *dialog,
                                        gint             response_id,
                                        gpointer         user_data)
{
 if(response_id == GTK_RESPONSE_CLOSE || 
    response_id == GTK_RESPONSE_CANCEL){
  gtk_widget_hide(aboutDialog);
 }
}

