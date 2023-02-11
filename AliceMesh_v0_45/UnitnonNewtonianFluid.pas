unit UnitnonNewtonianFluid;
// ћодуль дл€ задани€ параметров неньютоновских жидкостей.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, ExtCtrls, StdCtrls, ImgList, jpeg;

type
  TFormnonNewtonFluid = class(TForm)
    GBMain: TGroupBox;
    GBsellaw: TGroupBox;
    CBlaw: TComboBox;
    ImageOstwald_de_Vel: TImage;
    Imagecaisson: TImage;
    ImagePrandtl: TImage;
    ImageCarreau: TImage;
    ImagePowell_Eyring: TImage;
    ImageWilliamson: TImage;
    GBViscosityLimiter: TGroupBox;
    LabelMin: TLabel;
    Labelmax: TLabel;
    Editmin: TEdit;
    Editmax: TEdit;
    GBparam: TGroupBox;
    LabelA: TLabel;
    EditA: TEdit;
    LabelB: TLabel;
    EditB: TEdit;
    LabelC: TLabel;
    EditC: TEdit;
    ButtonApply: TButton;
    Labeln: TLabel;
    Editn: TEdit;
    procedure CBlawChange(Sender: TObject);
    procedure ButtonApplyClick(Sender: TObject);
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormnonNewtonFluid: TFormnonNewtonFluid;

implementation
 uses
      UnitUserDefinedFluidMaterial, VisualUnit;
{$R *.dfm}

procedure TFormnonNewtonFluid.CBlawChange(Sender: TObject);
begin
   // изменение закона приводит к изменению показа картинки
   // на которой написан закон

   ImageOstwald_de_Vel.Visible:=false;
   Imagecaisson.Visible:=false;
   ImagePrandtl.Visible:=false;
   ImageCarreau.Visible:=false;
   ImagePowell_Eyring.Visible:=false;
   ImageWilliamson.Visible:=false;

   case CBlaw.ItemIndex of
     0 : begin
            // power-law fluid
            // ќствальд де ¬ель
            ImageOstwald_de_Vel.Visible:=true;
            LabelB.Visible:=false;
            EditB.Visible:=false;
            LabelC.Visible:=false;
            EditC.Visible:=false;
            Labeln.Visible:=true;
            Editn.Visible:=true;
         end;
     1 : begin
            //  ессон
            Imagecaisson.Visible:=true;
            LabelB.Visible:=true;
            EditB.Visible:=true;
            LabelC.Visible:=false;
            EditC.Visible:=false;
            Labeln.Visible:=false;
            Editn.Visible:=false;
         end;
     2 : begin
            // ѕрандтль
            ImagePrandtl.Visible:=true;
            LabelB.Visible:=true;
            EditB.Visible:=true;
            LabelC.Visible:=false;
            EditC.Visible:=false;
            Labeln.Visible:=false;
            Editn.Visible:=false;
         end;
     3 : begin
            // Carreau
            ImageCarreau.Visible:=true;
            LabelB.Visible:=true;
            EditB.Visible:=true;
            LabelC.Visible:=true;
            EditC.Visible:=true;
            Labeln.Visible:=true;
            Editn.Visible:=true;
         end;
     4 : begin
            // ѕауэлл-Ёйринг
            ImagePowell_Eyring.Visible:=true;
            LabelB.Visible:=true;
            EditB.Visible:=true;
            LabelC.Visible:=true;
            EditC.Visible:=true;
            Labeln.Visible:=false;
            Editn.Visible:=false;
         end;
     5 : begin
            // ”иль€мсон
            ImageWilliamson.Visible:=true;
            LabelB.Visible:=true;
            EditB.Visible:=true;
            LabelC.Visible:=true;
            EditC.Visible:=true;
            Labeln.Visible:=false;
            Editn.Visible:=false;
         end;
   end;
end;

// считывание параметров модели неньютоновской жидкости.
procedure TFormnonNewtonFluid.ButtonApplyClick(Sender: TObject);
begin
   // —читывание закона изменени€ динамической в€зкости дл€ неньютоновских жидкостей
   Laplas.workmat[Laplas.itek].ilawmu:=CBlaw.ItemIndex+1; // номер закона дл€ динамической в€зкости
   Laplas.workmat[Laplas.itek].mumin:=StrToFloat(Editmin.Text); // минимальное значение динамической в€зкости
   Laplas.workmat[Laplas.itek].mumax:=StrToFloat(Editmax.Text); // максимальное значение динамической в€зкости
   Laplas.workmat[Laplas.itek].mu:=0.5*(Laplas.workmat[Laplas.itek].mumin+Laplas.workmat[Laplas.itek].mumax);
   Laplas.workmat[Laplas.itek].Amu:=StrToFloat(EditA.Text); // параметры
   Laplas.workmat[Laplas.itek].Bmu:=StrToFloat(EditB.Text); // моделей
   Laplas.workmat[Laplas.itek].Cmu:=StrToFloat(EditC.Text);
   Laplas.workmat[Laplas.itek].degreennmu:=StrToFloat(Editn.Text); // показатель степени дл€ некоторых моделей
end;  // non Newtonian Fluid

end.
