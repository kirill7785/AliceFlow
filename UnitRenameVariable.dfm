object FormRenameVar: TFormRenameVar
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Rename variable'
  ClientHeight = 153
  ClientWidth = 417
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 0
    Top = 0
    Width = 417
    Height = 153
    Caption = 'Select variable for rename'
    Color = clMoneyGreen
    ParentBackground = False
    ParentColor = False
    TabOrder = 0
    object LabelRenamecandidate: TLabel
      Left = 16
      Top = 27
      Width = 89
      Height = 13
      Caption = 'Rename candidate'
    end
    object Label1: TLabel
      Left = 24
      Top = 64
      Width = 290
      Height = 13
      Caption = 'Please, enter new name vriable. First simbol must be equal $'
    end
    object ComboBox1: TComboBox
      Left = 111
      Top = 24
      Width = 290
      Height = 21
      TabOrder = 0
      Text = 'ComboBox1'
    end
    object EditNewName: TEdit
      Left = 232
      Top = 80
      Width = 169
      Height = 21
      TabOrder = 1
    end
    object Button1: TButton
      Left = 296
      Top = 112
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 2
      OnClick = Button1Click
    end
  end
end
