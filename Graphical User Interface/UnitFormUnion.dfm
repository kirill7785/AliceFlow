object FormUnion: TFormUnion
  Left = 527
  Top = 153
  AutoSize = True
  Caption = 'Union'
  ClientHeight = 393
  ClientWidth = 305
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object PanelMain: TPanel
    Left = 0
    Top = 0
    Width = 305
    Height = 393
    Color = clMoneyGreen
    TabOrder = 0
    object Labelname: TLabel
      Left = 16
      Top = 16
      Width = 26
      Height = 13
      Caption = 'name'
    end
    object Labelmassa: TLabel
      Left = 120
      Top = 40
      Width = 56
      Height = 13
      Caption = 'Labelmassa'
    end
    object ButtonApply: TButton
      Left = 160
      Top = 360
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 0
      OnClick = ButtonApplyClick
    end
    object CheckBoxVisible: TCheckBox
      Left = 16
      Top = 40
      Width = 97
      Height = 17
      Caption = 'Visible Objects'
      Checked = True
      State = cbChecked
      TabOrder = 1
      OnClick = CheckBoxVisibleClick
    end
    object Editname: TEdit
      Left = 56
      Top = 8
      Width = 121
      Height = 21
      TabOrder = 2
    end
    object CBMeshassemblesSeparat: TCheckBox
      Left = 16
      Top = 72
      Width = 153
      Height = 17
      Caption = 'Mesh asembles separately'
      TabOrder = 3
      OnClick = CBMeshassemblesSeparatClick
    end
    object GBmeshasembles: TGroupBox
      Left = 8
      Top = 104
      Width = 281
      Height = 249
      Caption = 'Mesh asembles separately parameters'
      TabOrder = 4
      object GBMesh: TGroupBox
        Left = 16
        Top = 136
        Width = 129
        Height = 105
        Caption = 'Mesh count'
        TabOrder = 0
        object Linx: TLabel
          Left = 16
          Top = 24
          Width = 13
          Height = 13
          Caption = 'inx'
        end
        object Liny: TLabel
          Left = 16
          Top = 48
          Width = 13
          Height = 13
          Caption = 'iny'
        end
        object Linz: TLabel
          Left = 16
          Top = 72
          Width = 13
          Height = 13
          Caption = 'inz'
        end
        object CBinx: TComboBox
          Left = 40
          Top = 16
          Width = 73
          Height = 21
          TabOrder = 0
          Text = '23'
          Items.Strings = (
            '4'
            '5'
            '6'
            '7'
            '8'
            '9'
            '10'
            '11'
            '12'
            '13'
            '14'
            '15'
            '16'
            '17'
            '18'
            '19'
            '20'
            '21'
            '22'
            '23'
            '24'
            '25'
            '26'
            '27'
            '28'
            '29'
            '30'
            '31'
            '32'
            '33')
        end
        object CBiny: TComboBox
          Left = 40
          Top = 40
          Width = 73
          Height = 21
          ItemIndex = 19
          TabOrder = 1
          Text = '23'
          Items.Strings = (
            '4'
            '5'
            '6'
            '7'
            '8'
            '9'
            '10'
            '11'
            '12'
            '13'
            '14'
            '15'
            '16'
            '17'
            '18'
            '19'
            '20'
            '21'
            '22'
            '23'
            '24'
            '25'
            '26'
            '27'
            '28'
            '29'
            '30'
            '31'
            '32'
            '33')
        end
        object CBinz: TComboBox
          Left = 40
          Top = 64
          Width = 73
          Height = 21
          ItemIndex = 19
          TabOrder = 2
          Text = '23'
          Items.Strings = (
            '4'
            '5'
            '6'
            '7'
            '8'
            '9'
            '10'
            '11'
            '12'
            '13'
            '14'
            '15'
            '16'
            '17'
            '18'
            '19'
            '20'
            '21'
            '22'
            '23'
            '24'
            '25'
            '26'
            '27'
            '28'
            '29'
            '30'
            '31'
            '32'
            '33')
        end
      end
      object GBadditsize: TGroupBox
        Left = 16
        Top = 24
        Width = 249
        Height = 105
        Caption = 'Additional size'
        TabOrder = 1
        object LXmin: TLabel
          Left = 16
          Top = 24
          Width = 23
          Height = 13
          Caption = 'Xmin'
        end
        object LYmin: TLabel
          Left = 16
          Top = 48
          Width = 23
          Height = 13
          Caption = 'Ymin'
        end
        object LZmin: TLabel
          Left = 16
          Top = 72
          Width = 23
          Height = 13
          Caption = 'Zmin'
        end
        object LXmax: TLabel
          Left = 136
          Top = 24
          Width = 26
          Height = 13
          Caption = 'Xmax'
        end
        object LYmax: TLabel
          Left = 136
          Top = 48
          Width = 26
          Height = 13
          Caption = 'Ymax'
        end
        object LZmax: TLabel
          Left = 136
          Top = 72
          Width = 26
          Height = 13
          Caption = 'Zmax'
        end
        object Exmin: TEdit
          Left = 48
          Top = 16
          Width = 73
          Height = 21
          TabOrder = 0
        end
        object Eymin: TEdit
          Left = 48
          Top = 40
          Width = 73
          Height = 21
          TabOrder = 1
        end
        object Ezmin: TEdit
          Left = 48
          Top = 64
          Width = 73
          Height = 21
          TabOrder = 2
        end
        object EXmax: TEdit
          Left = 168
          Top = 16
          Width = 73
          Height = 21
          TabOrder = 3
        end
        object EYmax: TEdit
          Left = 168
          Top = 40
          Width = 73
          Height = 21
          TabOrder = 4
        end
        object EZmax: TEdit
          Left = 168
          Top = 64
          Width = 73
          Height = 21
          TabOrder = 5
        end
      end
    end
  end
end
