// Library under test.
//#include "/usr/local/include/gtkimageview/gtkanimview.h"
//#include "/usr/local/include/gtkimageview/gtkimagescrollwin.h"
//#include "/usr/local/include/gtkimageview/gtkimagetooldragger.h"
//#include "/usr/local/include/gtkimageview/gtkimagetoolpainter.h"
//#include "/usr/local/include/gtkimageview/gtkimagetoolselector.h"

#include <gtkimageview/gtkanimview.h>
#include <gtkimageview/gtkimagescrollwin.h>
#include <gtkimageview/gtkimagetooldragger.h>
#include <gtkimageview/gtkimagetoolpainter.h>
#include <gtkimageview/gtkimagetoolselector.h>

//#include <gtkanimview.h>
//#include <gtkimagescrollwin.h>
//#include <gtkimagetooldragger.h>
//#include <gtkimagetoolpainter.h>
//#include <gtkimagetoolselector.h>


// Defines for backwards compatibility with GTK+ 2.6
#if !GTK_CHECK_VERSION(2,8,0)
#define GTK_STOCK_FULLSCREEN ""
#endif

#define MULTIPOS_RUN     1
#define POS_RUN      2
#define SET_RUN      3

#define BF_FL_MAPPING_TIMETAG  0
#define BF_FL_MAPPING_FILENAME 1
#define BF_FL_MAPPING_LIST     2

#define BRIGHTFIELD             0
#define CONFOCAL_TRANSMISSION   1
#define METAMORPH_DECONVOLUTION 2
#define MEMBRANE_TAGGED_FL      3
#define UNKNOWN_IMAGE_TYPE      4

#define ALIGN_NO      0
#define ALIGN_FL_t0   1
#define ALIGN_FL_BF   2 

#define CELL_ALIGN_NO      0
#define CELL_ALIGN_CELL    1
#define CELL_ALIGN_BOUND   2 

#define FRET_NUCLEUS_TOP    0
#define FRET_NUCLEUS_BOTTOM 1



#define LIST_FILE 0  //List of phase and fuorescence images passed as arguments
                     //default in version 1.1
#define OIF_FILE 1   //Olymus Image File
#define OEM_FILE 2   //Open Enviroment for Microscopy file


typedef struct _params{
  //segmentation parameters //
  gdouble max_dist_over_waist,
          max_split_over_minor;

  gdouble background_reject_factor;
  gdouble tracking_comparison;  
  gint    min_pixels_per_cell,
          max_pixels_per_cell;

  gint       align_fl;     //align FL to X
  gint       align_cells;  //wiggling cells

  //only BRIGHTFIELD allowed in v0.1
  gint image_type; 

  gint splitted_fret_image,
       fret_nucleus;

  gint  bf_fl_mapping;

  gint  treat_brightfield_as_fluorescence;
 
  GString *output_basename,
          *output_path, //not settable. automatic
          *exec_file;
  
  GString *bf_token, *bf_path,
          *fl_token, *fl_path,
          *flat_token, *flat_path,
          *dark_token, *dark_path;
  GString *nucleus_channel; //third_image;
  gboolean bdark, bflat,bnucleus;
  gboolean bpos,btime,bsep;
  GString *pos_token, *time_token,*separator;        
  gboolean same_path;
  gboolean bPAW; 

  GString *tmp_path;

  //the following are not used in v0.1
  gdouble max_dist_over_waist_t0,
          max_split_over_minor_t0;
  gboolean   segment_values_t0;        
  gboolean   force_nucleus_in_center;
  GString   *mappingFilePath;
  gboolean  do_recomb;
  gboolean  all_recomb,
            nuc_recomb,
            fl_recomb;
  GString   *nuc_recomb_ch,
            *fl_recomb_ch;
  gdouble   nuc_recomb_val,
            fl_recomb_val;
  GString   *param_file;


}PARAMS;

typedef struct _img{
 gchar *fname,*fname_utf;
 gchar *pos,*pos_utf;
 gchar *time,*time_utf;
 gdouble timetag; //metamorph  diasJulianos + fraccion del dia en segundos.
}IMG;

typedef struct _imgList{
 GList   *gl;
 GString *token;
 GString *path;
}ImgList;

////////////////////////////////////////////
void InitParams(void);


void UpdateEntryFileName(void);
void Param2SegmentationDialog(void);
void Param2InputDialog(void);
void Param2ImageSetupDialog(void);
void Param2SetupDialog(void);

void InputDialog2Param(void);
void SegmentationDialog2Param(void);
void ImageSetupDialog2Param(void);
void SetupDialog2Param(void);

gboolean CheckSegmentationForm(void);
gboolean CheckInputForm(void);
gboolean CheckImageSetupForm(void);
gboolean CheckExecFile(void);

int  LoadParams(gchar *fname);
gint StatusMessage(gint cid,gchar *msg);

gboolean SetComboBoxItem(GtkComboBox *w,gchar *saux);
int      GetComboBoxItem(GtkComboBox *w,gchar **);
GtkResponseType message_dialog( GtkMessageType type,
                                GtkWindow *parent,
                                const gchar *title,
                                const gchar *message);
void ReportWarning(gchar *mesg);


gboolean LoadImageNames(void);
GList    *GetFileNames(gchar *ppath, gchar *ttoken);
GString  *GetFileNamePattern(gchar *imgType);


void CreateAndFillModel2(GList *bf,gboolean b);
void CreateAndFillModel1(GList *bf,gboolean b);
void CreateAndFillModelX(GList *bf, gchar *ftoken);
gboolean GetTreeIter(gchar *saux, 
                     GtkTreeIter *parent, 
                     GtkTreeIter *iter);
gint GetTreeDepth(GtkTreeView *tv);

gint GetClosestTimeTreeIter(gchar *pchar,
                                GtkTreeIter *parent,
                                GtkTreeIter *iter);
gint GetClosestTimeTagTreeIter(gdouble i ,
                                GtkTreeIter *parent,
                                GtkTreeIter *iter);                               

void LoadImages(gchar *name);

void  view_popup_menu (GtkWidget *treeview, GdkEventButton *event, gpointer userdata);
void  view_popup_menu_onDoTest (GtkWidget *menuitem, gpointer userdata);
void  view_popup_menu_onDoRun (GtkWidget *menuitem, gpointer userdata);
gboolean view_onPopupMenu (GtkWidget *treeview, gpointer userdata);
gboolean view_onButtonPressed (GtkWidget *treeview, GdkEventButton *event, gpointer userdata);

void TestCellID();
gboolean LaunchCellID(GtkTreeIter *p,gboolean btest,gint itype);
void WriteTXTfiles_Set(GtkTreeIter *pIterSet,gchar *wmode,gboolean btest);
int  SaveParams(gboolean btest);
gboolean filter_visible_function (GtkTreeModel *model, GtkTreeIter *iter, gpointer data);
gchar *InsertTestBFNode(GtkTreeIter *piter);

gint  sort_func_filename (GtkTreeModel *model,
                 GtkTreeIter  *a,
                 GtkTreeIter  *b,
                 gpointer      userdata);
gboolean GetBFTreeIter(GtkTreeIter *piter, GtkTreeIter *pIterSet);
void RunCellID(void);
void CleanTMPdir(void);
